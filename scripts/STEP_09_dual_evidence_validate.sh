#!/usr/bin/env bash
#SBATCH --job-name=bp_validate
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/bp_validate_%j.out
#SBATCH --error=logs/bp_validate_%j.err
# ============================================================================
# STEP_09_dual_evidence_validate.sh
#
# Reciprocal validation of breakpoint regions using two independent evidence
# types:
#
#   1. K-MER EVIDENCE  : mashmap multi-pi sweep (pi = 70, 80, 85, 90, 95)
#                        Each pi threshold gives blocks. Blocks at pi=95 are
#                        ironclad; blocks only visible at pi=70 are ancient.
#
#   2. PROTEIN EVIDENCE: miniprot alignment of a reference proteome against
#                        the same windows. Gene hits are sparse but robust to
#                        intron drift, synonymous substitutions, and repeats.
#
# Logic:
#   - If both evidences agree on orthology for a tile  -> HIGH_CONFIDENCE
#   - If only one evidence supports orthology          -> MEDIUM_CONFIDENCE
#   - If neither (or they disagree)                    -> LOW_CONFIDENCE
#
# Breakpoints with HIGH-confidence flanks on both sides are solid calls.
# Breakpoints with LOW-confidence flanks are assembly-artifact suspects or
# regions of extreme divergence worth investigating separately.
#
# Usage:
#   sbatch STEP_09_dual_evidence_validate.sh \
#       <bp_bed> <concat_fa> <reference_proteome.faa>
#
# The reference proteome should be a trusted gene set — e.g., your existing
# Cgar or Cmac protein predictions, or zebrafish RefSeq proteins for a more
# conserved anchor.
# ============================================================================
set -euo pipefail

BP_BED="${1:?Usage: STEP_09 <bp_bed> <concat_fa> <proteome.faa>}"
CONCAT_FA="${2:?Usage: STEP_09 <bp_bed> <concat_fa> <proteome.faa>}"
PROTEOME="${3:?Usage: STEP_09 <bp_bed> <concat_fa> <proteome.faa>}"

module load Miniconda3 || true
source activate assembly

for tool in samtools mashmap minimap2 miniprot; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found"; exit 1; }
done

OUTDIR="results/09_dual_evidence"
mkdir -p "$OUTDIR"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

# PI thresholds for the sweep
PI_LEVELS=(70 80 85 90 95)

# ============================================================================
# Phase 1: extract all breakpoint windows into one concatenated fasta
# ============================================================================
echo "[STEP_09] Phase 1: extract breakpoint windows"
WINDOWS_FA="$OUTDIR/all_windows.fa"
> "$WINDOWS_FA"
while IFS=$'\t' read -r region start end bp_id; do
    samtools faidx "$CONCAT_FA" "${region}:${start}-${end}" \
        | sed "s/^>.*/>${bp_id}__${region}_${start}_${end}/" \
        >> "$WINDOWS_FA"
done < "$BP_BED"
samtools faidx "$WINDOWS_FA"
echo "  $(grep -c '^>' "$WINDOWS_FA") windows extracted"

# ============================================================================
# Phase 2: mashmap multi-pi sweep (windows vs windows, all-vs-all)
# ============================================================================
echo ""
echo "[STEP_09] Phase 2: mashmap multi-pi sweep"
MASH_DIR="$OUTDIR/mashmap_sweep"
mkdir -p "$MASH_DIR"

for pi in "${PI_LEVELS[@]}"; do
    echo "  pi=$pi"
    mashmap \
        -r "$WINDOWS_FA" \
        -q "$WINDOWS_FA" \
        -s 5000 \
        --pi "$pi" \
        -f one-to-one \
        --kmerComplexity 0.5 \
        -t "$THREADS" \
        -o "$MASH_DIR/windows_pi${pi}.paf" \
        2> "$MASH_DIR/windows_pi${pi}.log"
done

# ============================================================================
# Phase 3: miniprot protein alignment
# ============================================================================
echo ""
echo "[STEP_09] Phase 3: miniprot protein-to-genome alignment"
MP_DIR="$OUTDIR/miniprot"
mkdir -p "$MP_DIR"

# Index the windows FASTA
miniprot -t "$THREADS" -d "$MP_DIR/windows.mpi" "$WINDOWS_FA" 2> "$MP_DIR/index.log"

# Align proteome; --gff for easy parsing
miniprot -Iut "$THREADS" --gff \
    "$MP_DIR/windows.mpi" "$PROTEOME" \
    > "$MP_DIR/proteome_vs_windows.gff" \
    2> "$MP_DIR/align.log"

echo "  $(grep -c '^##PAF' "$MP_DIR/proteome_vs_windows.gff" || echo 0) protein alignments"
echo "  $(awk '$3 == "mRNA"' "$MP_DIR/proteome_vs_windows.gff" | wc -l) mRNAs"

# ============================================================================
# Phase 4: tile each window into 5 kb bins and compute dual-evidence score
# ============================================================================
echo ""
echo "[STEP_09] Phase 4: dual-evidence scoring"
python3 "$(dirname "$0")/STEP_09b_score_dual_evidence.py" \
    --bp-bed "$BP_BED" \
    --mash-dir "$MASH_DIR" \
    --miniprot-gff "$MP_DIR/proteome_vs_windows.gff" \
    --pi-levels "${PI_LEVELS[@]}" \
    --out "$OUTDIR/" \
    --tile-bp 5000

echo ""
echo "[STEP_09] Done. Outputs:"
echo "  $OUTDIR/tile_evidence.tsv              — per-tile evidence scores"
echo "  $OUTDIR/breakpoint_confidence.tsv      — per-breakpoint confidence calls"
echo "  $OUTDIR/evidence_tracks.pdf            — per-breakpoint evidence plot"
