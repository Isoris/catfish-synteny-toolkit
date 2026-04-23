#!/usr/bin/env bash
#SBATCH --job-name=prep_genome
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-10
#SBATCH --output=logs/prep_genome_%A_%a.out
#SBATCH --error=logs/prep_genome_%A_%a.err
# ============================================================================
# STEP_00b_prepare_genomes.sh
#
# Per-species genome preparation using the EXISTING ScaffoldKit Module 1
# helpers from the MODULE_CONSERVATION project. Reuses tested code rather
# than reimplementing.
#
# For each species in the manifest, produces:
#   <species>.scaffolds.fa           Picard-normalized, cleaned headers
#   <species>.scaffolds.fa.fai       samtools index
#   <species>.scaffolds.sizes        seqkit fx2tab
#   <species>.scaffolds.gaps         detgaps
#   <species>.scaffolds.ragtag.agp   RagTag splitasm AGP
#   <species>.contigs.fa             contigs after gap splitting
#   <species>.scaffolds.telomeres.tsv   Quartet TeloExplorer
#   <species>.contigs.telomeres.tsv
#   <species>.report.txt             human-readable summary
#   <species>.md5                    checksums
#
# These files are then consumed by STEP_00_prepare_panSN.sh which takes the
# .scaffolds.fa as input and concatenates with PanSN naming.
#
# Usage:
#   Set SCAFFOLDKIT_HELPERS env var to point to your scaffoldkit_module1 dir,
#   or edit the default below.
#
#   sbatch --array=0-$((N_SPECIES-1)) STEP_00b_prepare_genomes.sh [tier]
#
# If invoked without SLURM_ARRAY_TASK_ID (i.e. interactively), processes ALL
# species in the tier sequentially.
# ============================================================================
set -euo pipefail

TIER="${1:-all}"

module load Miniconda3 || true
source activate assembly

# ---- Locate scaffoldkit helpers ----
# Default path on LANTA from the MODULE_CONSERVATION project
: "${SCAFFOLDKIT_HELPERS:=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/MODULE_CONSERVATION/helpers/scaffoldkit_module1}"

if [[ ! -d "$SCAFFOLDKIT_HELPERS" ]]; then
    echo "ERROR: SCAFFOLDKIT_HELPERS not found at: $SCAFFOLDKIT_HELPERS"
    echo "Set the env variable: export SCAFFOLDKIT_HELPERS=/path/to/scaffoldkit_module1"
    exit 1
fi
echo "[STEP_00b] ScaffoldKit helpers: $SCAFFOLDKIT_HELPERS"

# ---- Tier filter ----
case "$TIER" in
    core)    TIERS_REGEX="^(core)$" ;;
    clarias) TIERS_REGEX="^(core|clarias_context)$" ;;
    all)     TIERS_REGEX="^(core|clarias_context|far_outgroup)$" ;;
    deep)    TIERS_REGEX="^(core|clarias_context|deep_outgroup)$" ;;
    *) echo "ERROR: tier must be core | clarias | all | deep"; exit 1 ;;
esac

MANIFEST="config/species_manifest.tsv"
OUTDIR_BASE="results/00_prepared_genomes"
LOGDIR="logs"
mkdir -p "$OUTDIR_BASE" "$LOGDIR"

# ---- Build list of species to process for this tier ----
mapfile -t SPECIES_LINES < <(tail -n +2 "$MANIFEST" | awk -F'\t' -v pat="$TIERS_REGEX" '$6 ~ pat && $1 != "" && $1 !~ /^#/')
N="${#SPECIES_LINES[@]}"
echo "[STEP_00b] Tier '$TIER': $N species to prepare"

# ---- Determine which indices to run ----
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # Running as array: process exactly one species
    INDICES=("$SLURM_ARRAY_TASK_ID")
    if (( SLURM_ARRAY_TASK_ID >= N )); then
        echo "  Array task $SLURM_ARRAY_TASK_ID >= $N species. Nothing to do."
        exit 0
    fi
else
    # Interactive: sweep all
    INDICES=()
    for ((i = 0; i < N; i++)); do INDICES+=("$i"); done
fi

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ============================================================================
# Per-species preparation function
# ============================================================================
prepare_one() {
    local sp_id="$1"
    local raw_fa="$2"
    local outdir="$3"

    mkdir -p "$outdir"
    local base="$sp_id"
    local scaff_fa="${outdir}/${base}.scaffolds.fa"

    # ── B1. Normalize FASTA (Picard) ──
    if [[ ! -f "${scaff_fa}" ]]; then
        echo "    [B1] Normalizing FASTA for $sp_id..."
        if command -v picard &>/dev/null && [[ -f "${SCAFFOLDKIT_HELPERS}/module-1-datasets_normalize_fasta.sh" ]]; then
            bash "${SCAFFOLDKIT_HELPERS}/module-1-datasets_normalize_fasta.sh" \
                --input "${raw_fa}" \
                --output "${scaff_fa}" \
                --fasta-line-length 80 \
                --truncate-header-whitespace \
                >> "${LOGDIR}/00b_${sp_id}.log" 2>&1 || true
            # Header cleanup: strip trailing info, strip pipe prefixes
            sed -i -E 's/^>[^|]*\|/>/; s/^(>[^ ]+).*/\1/' "${scaff_fa}"
        elif command -v picard &>/dev/null; then
            echo "      Picard available but scaffoldkit helper not found, using Picard direct"
            picard NormalizeFasta \
                -I "${raw_fa}" -O "${scaff_fa}" \
                -LINE_LENGTH 80 \
                -TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE true \
                -QUIET true -VERBOSITY ERROR
            sed -i -E 's/^>[^|]*\|/>/; s/^(>[^ ]+).*/\1/' "${scaff_fa}"
        else
            echo "      [WARN] Picard not available — using awk fallback"
            awk '/^>/{print $1; next}{print}' "${raw_fa}" > "${scaff_fa}"
        fi
        samtools faidx "${scaff_fa}"
    fi

    # ── B2. Scaffold sizes (seqkit fx2tab) ──
    local sizes="${outdir}/${base}.scaffolds.sizes"
    if [[ ! -f "$sizes" ]]; then
        echo "    [B2] Scaffold sizes..."
        if command -v seqkit &>/dev/null; then
            seqkit fx2tab -nl "${scaff_fa}" > "${sizes}"
        else
            awk '{print $1"\t"$2}' "${scaff_fa}.fai" > "${sizes}"
        fi
    fi

    # ── B3. Gap detection (detgaps or Python fallback) ──
    local gaps="${outdir}/${base}.scaffolds.gaps"
    if [[ ! -f "$gaps" ]]; then
        echo "    [B3] Detecting gaps..."
        if command -v detgaps &>/dev/null; then
            detgaps "${scaff_fa}" > "${gaps}"
        else
            python3 - "${scaff_fa}" > "${gaps}" <<'PYEOF'
import sys
with open(sys.argv[1]) as f:
    name, pos, in_gap, gap_start = "", 0, False, 0
    for line in f:
        if line.startswith(">"):
            if in_gap: print(f"{name}\t{gap_start}\t{pos}")
            name = line[1:].split()[0]; pos = 0; in_gap = False
        else:
            for c in line.strip():
                if c in "Nn":
                    if not in_gap: gap_start = pos; in_gap = True
                else:
                    if in_gap: print(f"{name}\t{gap_start}\t{pos}"); in_gap = False
                pos += 1
    if in_gap: print(f"{name}\t{gap_start}\t{pos}")
PYEOF
        fi
    fi

    # ── B4. AGP via RagTag splitasm (optional — only if ragtag available) ──
    local agp="${outdir}/${base}.scaffolds.ragtag.agp"
    local contigs_fa="${outdir}/${base}.contigs.fa"
    if [[ ! -f "$agp" ]] && command -v ragtag.py &>/dev/null; then
        echo "    [B4] RagTag splitasm for AGP + contigs..."
        ragtag.py splitasm "${scaff_fa}" -o "${outdir}/ragtag_tmp" >> "${LOGDIR}/00b_${sp_id}.log" 2>&1 || true
        [[ -f "${outdir}/ragtag_tmp/split.agp" ]] && mv "${outdir}/ragtag_tmp/split.agp" "${agp}"
        [[ -f "${outdir}/ragtag_tmp/split.fa" ]]  && mv "${outdir}/ragtag_tmp/split.fa"  "${contigs_fa}"
        rm -rf "${outdir}/ragtag_tmp"
    fi

    # ── B5. Telomere detection (optional — only if Quartet available) ──
    local telos="${outdir}/${base}.scaffolds.telomeres.tsv"
    if [[ ! -f "$telos" ]] && command -v quartet_teloexplorer.py &>/dev/null; then
        echo "    [B5] Telomere analysis..."
        local cwd="$PWD"
        cd "$outdir"
        quartet_teloexplorer.py -i "$(basename "${scaff_fa}")" -c animal -m 100 \
            --noplot -p "${base}.scaffolds" \
            >> "${cwd}/${LOGDIR}/00b_${sp_id}.log" 2>&1 || true
        if [[ -f "${base}.scaffolds.telo.info" && -f "${SCAFFOLDKIT_HELPERS}/SCAFFOLDKIT_process_telos.sh" ]]; then
            bash "${SCAFFOLDKIT_HELPERS}/SCAFFOLDKIT_process_telos.sh" \
                "$(basename "${scaff_fa}")" "${base}.scaffolds" \
                >> "${cwd}/${LOGDIR}/00b_${sp_id}.log" 2>&1 || true
            [[ -f "${base}.scaffolds.telo.cleaned.info.tsv" ]] && \
                mv "${base}.scaffolds.telo.cleaned.info.tsv" "${telos}"
        fi
        # tidy temp files
        mkdir -p extra
        mv ${base}.scaffolds*telo* extra/ 2>/dev/null || true
        mv ${base}.scaffolds*_telos.tsv extra/ 2>/dev/null || true
        cd "$cwd"
    fi

    # ── B6. Seqkit filter contigs < MIN_CONTIG_BP (user mentioned 5 Mb filter) ──
    #     Only filter if the assembly has many small unplaced contigs; harmless otherwise.
    local filtered_fa="${outdir}/${base}.scaffolds.filtered.fa"
    if [[ ! -f "$filtered_fa" ]] && command -v seqkit &>/dev/null; then
        echo "    [B6] Filtering contigs < 5 Mb..."
        seqkit seq -m 5000000 "${scaff_fa}" > "${filtered_fa}" 2>/dev/null || true
        if [[ -s "$filtered_fa" ]]; then
            samtools faidx "${filtered_fa}"
        else
            echo "      WARNING: filtered FASTA empty — keeping all contigs"
            rm -f "${filtered_fa}"
            cp "${scaff_fa}" "${filtered_fa}"
            samtools faidx "${filtered_fa}"
        fi
    fi

    # ── B7. Report ──
    local report="${outdir}/${base}.report.txt"
    {
        echo "=========================================="
        echo "GENOME PREP REPORT — $sp_id"
        echo "=========================================="
        echo "Input:     $raw_fa"
        echo "Scaffolds: ${scaff_fa}"
        echo ""
        if [[ -f "${scaff_fa}.fai" ]]; then
            local total_bp n_seq largest
            total_bp=$(awk '{s+=$2}END{print s}' "${scaff_fa}.fai")
            n_seq=$(wc -l < "${scaff_fa}.fai")
            largest=$(sort -k2 -rn "${scaff_fa}.fai" | head -1 | awk '{print $1" ("$2" bp)"}')
            echo "Total bp:  $total_bp"
            echo "Sequences: $n_seq"
            echo "Largest:   $largest"
        fi
        [[ -f "$gaps" ]]  && echo "Gaps:      $(wc -l < "$gaps") regions"
        [[ -f "$agp" ]]   && echo "AGP:       $agp"
        [[ -f "$telos" ]] && echo "Telomeres: $telos"
        if [[ -f "$filtered_fa.fai" ]]; then
            echo ""
            echo "After <5Mb filter:"
            local f_total f_n
            f_total=$(awk '{s+=$2}END{print s}' "${filtered_fa}.fai")
            f_n=$(wc -l < "${filtered_fa}.fai")
            echo "  Total bp:  $f_total"
            echo "  Sequences: $f_n"
        fi
    } > "${report}"

    # ── B8. MD5 ──
    local md5="${outdir}/${base}.md5"
    md5sum "${scaff_fa}" "${filtered_fa}" 2>/dev/null > "${md5}" || true

    echo "  [DONE] $sp_id"
    cat "${report}"
}

# ============================================================================
# Process selected species
# ============================================================================
for idx in "${INDICES[@]}"; do
    IFS=$'\t' read -r sp_id genus species prefix fasta tier notes <<< "${SPECIES_LINES[$idx]}"

    echo ""
    echo "============================================================"
    echo "Preparing [$idx] $sp_id ($genus $species, tier=$tier)"
    echo "============================================================"

    if [[ ! -f "$fasta" ]]; then
        echo "  ERROR: $fasta not found — skipping $sp_id"
        continue
    fi

    prepare_one "$sp_id" "$fasta" "${OUTDIR_BASE}/${sp_id}"
done

echo ""
echo "[STEP_00b] Genome preparation complete for tier=$TIER"
echo "[STEP_00b] Per-species outputs: ${OUTDIR_BASE}/<species>/"
echo ""
echo "Next step: update config/species_manifest.tsv 'fasta_path' column to"
echo "point to the filtered FASTAs (${OUTDIR_BASE}/<species>/<species>.scaffolds.filtered.fa),"
echo "then run STEP_00_prepare_panSN.sh."
