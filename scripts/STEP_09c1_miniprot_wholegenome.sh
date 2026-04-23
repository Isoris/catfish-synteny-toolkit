#!/usr/bin/env bash
#SBATCH --job-name=miniprot_wg
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --array=0-8       # adjust to n_species - 1
#SBATCH --output=logs/miniprot_wg_%A_%a.out
#SBATCH --error=logs/miniprot_wg_%A_%a.err
# ============================================================================
# STEP_09c1_miniprot_wholegenome.sh
#
# COLLECT phase for STEP_09c flank coherence analysis.
#
# Runs miniprot with LENIENT parameters against each species' full genome.
# Output: one GFF per species with all protein family hits.
#
# This is the "cast wide net" step — the collect-then-refine principle requires
# that collection is lenient enough to capture everything interesting.
# --outs=0.5 keeps hits at 50% of best score (default 0.95 is too strict for
# cross-species annotation transfer at ~90-95% ANI).
#
# SLURM array: one task per species in the manifest.
#
# Usage:
#   sbatch --array=0-$((N-1)) STEP_09c1_miniprot_wholegenome.sh <proteome.faa> <tier>
# ============================================================================
set -euo pipefail

PROTEOME="${1:?Usage: STEP_09c1 <proteome.faa> <tier>}"
TIER="${2:?Usage: STEP_09c1 <proteome.faa> <tier>}"

module load Miniconda3 || true
source activate assembly
command -v miniprot >/dev/null 2>&1 || { echo "ERROR: miniprot not found"; exit 1; }

OUTDIR="results/09c_wg_miniprot"
mkdir -p "$OUTDIR"
THREADS="${SLURM_CPUS_PER_TASK:-64}"

MANIFEST="config/species_manifest.tsv"

# Select species for this array task, filtered by tier
case "$TIER" in
    core)    TIERS_REGEX="^(core)$" ;;
    clarias) TIERS_REGEX="^(core|clarias_context)$" ;;
    all)     TIERS_REGEX="^(core|clarias_context|far_outgroup)$" ;;
    deep)    TIERS_REGEX="^(core|clarias_context|deep_outgroup)$" ;;
    *) echo "ERROR: tier must be core | clarias | all | deep"; exit 1 ;;
esac

# Build species list from manifest
mapfile -t SPECIES_LINES < <(tail -n +2 "$MANIFEST" | awk -F'\t' -v pat="$TIERS_REGEX" '$6 ~ pat && $1 != "" && $1 !~ /^#/')

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
if (( TASK_ID >= ${#SPECIES_LINES[@]} )); then
    echo "Task $TASK_ID out of range (max ${#SPECIES_LINES[@]} species)"
    exit 0
fi

IFS=$'\t' read -r sp_id genus species prefix fasta tier notes <<< "${SPECIES_LINES[$TASK_ID]}"

echo "[STEP_09c1] Task $TASK_ID: $sp_id ($genus $species, tier=$tier)"
echo "[STEP_09c1] Fasta:  $fasta"
echo "[STEP_09c1] Lenient collection (--outs=0.5), threads=$THREADS"

if [[ ! -f "$fasta" ]]; then
    echo "ERROR: $fasta not found"
    exit 1
fi

OUT_GFF="$OUTDIR/${sp_id}.miniprot.gff"
OUT_LOG="$OUTDIR/${sp_id}.miniprot.log"

# Lenient parameters:
#   -I        : auto intron size based on genome size
#   -u        : output unmapped queries too (for completeness stats)
#   --outs=0.5: keep suboptimal hits at >=50% of best score
#   --gff     : GFF3 output
time miniprot \
    -I -u -t "$THREADS" \
    --outs=0.5 \
    --gff \
    "$fasta" "$PROTEOME" \
    > "$OUT_GFF" 2> "$OUT_LOG"

n_mrna=$(awk '$3 == "mRNA"' "$OUT_GFF" | wc -l)
echo "[STEP_09c1] Done: $n_mrna mRNAs in $OUT_GFF"
