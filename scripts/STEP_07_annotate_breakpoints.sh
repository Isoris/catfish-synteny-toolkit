#!/usr/bin/env bash
#SBATCH --job-name=bp_annotate
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/bp_annotate_%j.out
#SBATCH --error=logs/bp_annotate_%j.err
# ============================================================================
# STEP_07_annotate_breakpoints.sh
#
# STAGE C — Local annotation at each breakpoint window.
#
# For each breakpoint interval from STEP_05 (breakpoints.bed), extract the
# sequence window and run a panel of established annotation tools:
#
#   TRF  (Benson 1999)            : tandem repeats + microsatellites
#   IRF  (Benson 2002)            : inverted repeats / dyad symmetry  <-- Kuang Fig 2C
#   minimap2 self-align           : segmental duplications within window
#   user --TE-lib (RepeatMasker)  : transposable element occupancy
#
# Output per-breakpoint:
#   bp_<id>/
#     window.fa                 : extracted sequence ±100 kb
#     trf.dat                   : tandem repeats
#     irf.dat                   : inverted repeats (palindromes)
#     selfaln.paf               : minimap2 self-alignment (SD detection)
#     te.out                    : RepeatMasker output if TE lib provided
#     annotations.bed           : merged BED of all features
#
# Plus a summary table across all breakpoints.
#
# Usage:
#   sbatch STEP_07_annotate_breakpoints.sh <bp_bed> <concat_fa> [TE_lib]
# ============================================================================
set -euo pipefail

BP_BED="${1:?Usage: STEP_07_annotate_breakpoints.sh <bp_bed> <concat_fa> [TE_lib]}"
CONCAT_FA="${2:?Usage: STEP_07_annotate_breakpoints.sh <bp_bed> <concat_fa> [TE_lib]}"
TE_LIB="${3:-}"

module load Miniconda3 || true
source activate assembly

# Required tools: samtools, trf, irf, minimap2, RepeatMasker (if TE lib given)
for tool in samtools trf minimap2; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found in PATH"; exit 1; }
done

if command -v irf >/dev/null 2>&1; then
    IRF_AVAILABLE=1
else
    echo "WARNING: irf not found — dyad symmetry analysis will be skipped."
    echo "         Install IRF from https://tandem.bu.edu/irf/irf.download.html"
    IRF_AVAILABLE=0
fi

if [[ -n "$TE_LIB" ]]; then
    command -v RepeatMasker >/dev/null 2>&1 || { echo "ERROR: RepeatMasker not found but TE_LIB provided"; exit 1; }
fi

OUTDIR="results/07_bp_annotation"
mkdir -p "$OUTDIR"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

# Summary table header
SUMMARY="$OUTDIR/breakpoint_annotation_summary.tsv"
printf "bp_id\tregion\tlength_bp\ttrf_n_repeats\ttrf_satellite_bp\ttrf_microsat_bp\tirf_n_palindromes\tirf_max_stem_bp\tselfaln_n_hits\tselfaln_max_len\tte_pct\n" > "$SUMMARY"

# ============================================================================
# Loop over breakpoints
# ============================================================================
n_total=$(wc -l < "$BP_BED")
i=0
while IFS=$'\t' read -r region start end bp_id; do
    i=$((i + 1))
    echo ""
    echo "============================================================"
    echo "[$i/$n_total] $bp_id : $region:$start-$end"
    echo "============================================================"

    bp_dir="$OUTDIR/$bp_id"
    mkdir -p "$bp_dir"
    window_fa="$bp_dir/window.fa"

    # ---- Extract window ----
    samtools faidx "$CONCAT_FA" "${region}:${start}-${end}" > "$window_fa"
    length=$((end - start))

    # ---- TRF: tandem repeats ----
    # Recommended params: 2 7 7 80 10 50 500
    #   Match=2, Mismatch=7, Delta=7, PM=80, PI=10, Minscore=50, MaxPeriod=500
    # -d = produce data file, -h = suppress html, -ngs = compact output
    echo "  [TRF] tandem repeats..."
    (cd "$bp_dir" && trf "$(basename "$window_fa")" 2 7 7 80 10 50 500 -d -h -ngs > trf.dat 2>trf.log) || true

    # Parse TRF output
    # ngs format: >seqname\n<start> <end> <period> <copies> <consensus_size> <pid> <indel%> <score> <A%> <C%> <G%> <T%> <entropy> <consensus> <sequence>
    trf_n=$(grep -v "^>" "$bp_dir/trf.dat" 2>/dev/null | awk 'NF>0' | wc -l)
    # Satellite bp = sum of lengths where period >= 100
    trf_sat_bp=$(grep -v "^>" "$bp_dir/trf.dat" 2>/dev/null | awk 'NF>0 && $3>=100 {s += $2-$1+1} END{print s+0}')
    # Microsat bp = sum of lengths where period 1-6
    trf_micro_bp=$(grep -v "^>" "$bp_dir/trf.dat" 2>/dev/null | awk 'NF>0 && $3>=1 && $3<=6 {s += $2-$1+1} END{print s+0}')

    # ---- IRF: inverted repeats (dyad symmetry, Kuang Fig 2C equivalent) ----
    irf_n=0
    irf_max_stem=0
    if [[ "$IRF_AVAILABLE" == "1" ]]; then
        echo "  [IRF] inverted repeats / palindromes..."
        # IRF params: Match Mismatch Delta PM PI Minscore Maxlength MaxLoop
        # Recommended: 2 3 5 80 10 40 500000 10000
        (cd "$bp_dir" && irf "$(basename "$window_fa")" 2 3 5 80 10 40 500000 10000 -d -h 2>irf.log) || true

        # IRF .dat file (find it — naming is complex)
        irf_dat=$(ls "$bp_dir"/*.irf.*.dat 2>/dev/null | head -1 || echo "")
        if [[ -n "$irf_dat" && -f "$irf_dat" ]]; then
            cp "$irf_dat" "$bp_dir/irf.dat"
            # Count palindromes (lines with content after header)
            irf_n=$(awk 'NF>=8 && $1 ~ /^[0-9]+$/' "$bp_dir/irf.dat" | wc -l)
            irf_max_stem=$(awk 'NF>=8 && $1 ~ /^[0-9]+$/ {len=$2-$1+1; if(len>max)max=len} END{print max+0}' "$bp_dir/irf.dat")
        fi
    fi

    # ---- minimap2 self-alignment: segmental duplications within window ----
    echo "  [minimap2] self-alignment for SDs..."
    minimap2 -X -cx asm20 -t "$THREADS" "$window_fa" "$window_fa" 2>/dev/null > "$bp_dir/selfaln.paf" || true
    # Count hits that are NOT identity-line (same coordinates) and >=1kb
    sd_n=$(awk 'BEGIN{FS="\t"} {
        if ($1 == $6 && $3 == $8 && $4 == $9) next;  # skip identity
        if ($11 >= 1000) print
    }' "$bp_dir/selfaln.paf" | wc -l)
    sd_maxlen=$(awk 'BEGIN{FS="\t"; m=0} {
        if ($1 == $6 && $3 == $8 && $4 == $9) next;
        if ($11 > m) m = $11
    } END{print m+0}' "$bp_dir/selfaln.paf")

    # ---- TE annotation (optional) ----
    te_pct=0
    if [[ -n "$TE_LIB" ]]; then
        echo "  [RepeatMasker] TE annotation..."
        (cd "$bp_dir" && RepeatMasker -lib "$TE_LIB" -pa "$THREADS" -nolow -no_is -q window.fa >rm.log 2>&1) || true
        if [[ -f "$bp_dir/window.fa.out" ]]; then
            # Sum masked bases from .out
            te_bp=$(awk 'NR>3 {s += $7-$6+1} END{print s+0}' "$bp_dir/window.fa.out")
            te_pct=$(awk -v t="$te_bp" -v l="$length" 'BEGIN{printf "%.2f", 100*t/l}')
        fi
    fi

    # ---- Write row to summary ----
    printf "%s\t%s:%d-%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" \
        "$bp_id" "$region" "$start" "$end" "$length" \
        "$trf_n" "$trf_sat_bp" "$trf_micro_bp" \
        "$irf_n" "$irf_max_stem" \
        "$sd_n" "$sd_maxlen" "$te_pct" >> "$SUMMARY"

    echo "  done: $trf_n tandem, $irf_n palindromes, $sd_n SDs, ${te_pct}% TE"

done < "$BP_BED"

echo ""
echo "[STEP_07] Done. Annotation summary: $SUMMARY"
