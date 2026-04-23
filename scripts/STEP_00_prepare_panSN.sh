#!/usr/bin/env bash
# ============================================================================
# STEP_00_prepare_panSN.sh
#   Reads config/species_manifest.tsv and concatenates selected genomes into
#   a single PanSN-named fasta for wfmash.
#
#   PanSN spec: >{prefix}#{haplotype}#{original_contig_name}
#   The '#' delimiter lets wfmash's -Y flag skip self-species mappings.
#
# Usage:
#   bash STEP_00_prepare_panSN.sh <tier_filter>
#   tier_filter: "core" | "clarias" | "all"
# ============================================================================
set -euo pipefail

TIER_FILTER="${1:-core}"
MANIFEST="config/species_manifest.tsv"
OUTDIR="results/00_panSN_inputs"
mkdir -p "$OUTDIR"

case "$TIER_FILTER" in
    core)     TIERS_REGEX="^(core)$" ;;
    clarias)  TIERS_REGEX="^(core|clarias_context)$" ;;
    all)      TIERS_REGEX="^(core|clarias_context|far_outgroup)$" ;;
    deep)     TIERS_REGEX="^(core|clarias_context|deep_outgroup)$" ;;
    *) echo "ERROR: tier_filter must be one of: core | clarias | all | deep"; exit 1 ;;
esac

OUT_FA="${OUTDIR}/catfish_${TIER_FILTER}.fa"
> "$OUT_FA"

echo "[STEP_00] Tier filter: $TIER_FILTER"
echo "[STEP_00] Output FASTA: $OUT_FA"
echo ""

n_species=0
while IFS=$'\t' read -r sp_id genus species prefix fasta tier notes; do
    [[ "$sp_id" =~ ^# ]] && continue  # skip comments
    [[ -z "$sp_id" ]] && continue
    [[ ! "$tier" =~ $TIERS_REGEX ]] && continue

    # Prefer the prepared (normalized + filtered) FASTA from STEP_00b if it exists
    prepared="results/00_prepared_genomes/${sp_id}/${sp_id}.scaffolds.filtered.fa"
    if [[ -f "$prepared" ]]; then
        echo "  [use prepared] $sp_id -> $prepared"
        fasta="$prepared"
    fi

    if [[ ! -f "$fasta" ]]; then
        echo "  WARNING: $fasta not found, skipping $sp_id"
        continue
    fi

    echo "  Adding $sp_id ($genus $species, tier=$tier)"

    # Rewrite headers to PanSN format: >prefix#1#origname
    if [[ "$fasta" == *.gz ]]; then
        zcat "$fasta"
    else
        cat "$fasta"
    fi | awk -v p="$prefix" '
        /^>/ {
            name = substr($0, 2)
            # strip anything after first whitespace
            sub(/[[:space:]].*$/, "", name)
            print ">" p "#1#" name
            next
        }
        { print }
    ' >> "$OUT_FA"
    n_species=$((n_species + 1))
done < <(tail -n +2 "$MANIFEST")

echo ""
echo "[STEP_00] Concatenated $n_species genomes"
echo "[STEP_00] Compressing with bgzip..."

# wfmash requires bgzip-compressed input with .fai index for efficiency
bgzip -f -@ 16 "$OUT_FA"
samtools faidx "${OUT_FA}.gz"

echo "[STEP_00] Done. Output files:"
ls -lh "${OUT_FA}.gz"*
