#!/usr/bin/env bash
#SBATCH --job-name=wfmash_catfish
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=logs/wfmash_%j.out
#SBATCH --error=logs/wfmash_%j.err
# ============================================================================
# STEP_01_wfmash_allvsall.sh
#   Run wfmash all-vs-all on PanSN-concatenated catfish genomes.
#   Produces:
#     - <tier>.paf         : base-accurate alignments (if not -m mode)
#     - <tier>.scaffolds.paf : filtered synteny scaffolds (the key output for breakpoints)
#
# Usage:
#   sbatch STEP_01_wfmash_allvsall.sh <tier> [mode]
#     tier: core | clarias | all
#     mode: approx (default, fast) | full (with WFA base alignment)
# ============================================================================
set -euo pipefail

TIER="${1:-core}"
MODE="${2:-approx}"
mkdir -p logs results/01_wfmash

module load Miniconda3 || true
source activate assembly

INPUT_FA="results/00_panSN_inputs/catfish_${TIER}.fa.gz"
OUT_PREFIX="results/01_wfmash/catfish_${TIER}"

if [[ ! -f "$INPUT_FA" ]]; then
    echo "ERROR: $INPUT_FA not found. Run STEP_00_prepare_panSN.sh first."
    exit 1
fi

# ----------------------------------------------------------------------------
# ANI threshold selection strategy:
#   - core tier (Clarias only):     ani50-3  (~90% ANI floor, Clarias-vs-Clarias)
#   - clarias tier:                 ani50-5  (slightly looser for Clarias context)
#   - all tier (with Silurus etc):  ani50-15 (must go lower to catch outgroup)
# ----------------------------------------------------------------------------
case "$TIER" in
    core)    ANI_PRESET="ani50-3"  ; SEG=50000  ; BLOCK=300000 ; SCAFF=100000 ;;
    clarias) ANI_PRESET="ani50-5"  ; SEG=50000  ; BLOCK=300000 ; SCAFF=100000 ;;
    all)     ANI_PRESET="ani50-15" ; SEG=30000  ; BLOCK=200000 ; SCAFF=80000  ;;
    deep)    ANI_PRESET="70"       ; SEG=20000  ; BLOCK=100000 ; SCAFF=50000  ;;
esac

# ----------------------------------------------------------------------------
# Common flags:
#   -Y '#'   : PanSN delimiter, skip intra-species mappings
#   -F 5e-5  : aggressive frequent-kmer filter (the flag you were remembering)
#   -n 1     : one mapping per segment (true synteny, not multi-mapping)
#   -o       : strict one-to-one plane-sweep filter
#   -s SEG   : segment length (Kuang et al. used effectively ~1 Mb at chr scale;
#              50 kb gives finer breakpoints without noise explosion)
#   -l BLOCK : minimum block length (~3x SEG, sharpens large syntenic regions)
#   -S SCAFF : scaffold mass — minimum chain length to form a synteny scaffold
#   -j 500k  : scaffold jump — max gap within a scaffold
#   --scaffold-out : THIS is the key output for breakpoint detection
# ----------------------------------------------------------------------------

APPROX_FLAG=""
OUT_PAF="${OUT_PREFIX}.paf"
if [[ "$MODE" == "approx" ]]; then
    APPROX_FLAG="-m"
    OUT_PAF="${OUT_PREFIX}.approx.paf"
fi

echo "=============================================================="
echo "wfmash all-vs-all catfish synteny"
echo "  Tier:       $TIER"
echo "  Mode:       $MODE"
echo "  Input:      $INPUT_FA"
echo "  ANI preset: $ANI_PRESET"
echo "  Segment:    $SEG bp"
echo "  Block min:  $BLOCK bp"
echo "  Scaffold:   $SCAFF bp mass"
echo "  Threads:    ${SLURM_CPUS_PER_TASK:-64}"
echo "=============================================================="
echo ""

time wfmash \
    -t "${SLURM_CPUS_PER_TASK:-64}" \
    -Y '#' \
    $APPROX_FLAG \
    -p "$ANI_PRESET" \
    -s "$SEG" \
    -l "$BLOCK" \
    -n 1 \
    -o \
    -F 0.00005 \
    -S "$SCAFF" \
    -j 500000 \
    --scaffold-out "${OUT_PREFIX}.scaffolds.paf" \
    "$INPUT_FA" \
    > "$OUT_PAF"

echo ""
echo "[STEP_01] Done."
echo "  Alignments:          $OUT_PAF ($(wc -l < "$OUT_PAF") records)"
echo "  Synteny scaffolds:   ${OUT_PREFIX}.scaffolds.paf ($(wc -l < "${OUT_PREFIX}.scaffolds.paf") records)"
echo ""
echo "Next: python3 scripts/STEP_02_call_breakpoints.py --paf ${OUT_PREFIX}.scaffolds.paf"
