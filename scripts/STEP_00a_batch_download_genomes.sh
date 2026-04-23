#!/usr/bin/env bash
#SBATCH --job-name=ncbi_download
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/ncbi_download_%j.out
#SBATCH --error=logs/ncbi_download_%j.err
# ============================================================================
# STEP_00a_batch_download_genomes.sh
#
# Reads a TSV of accessions and downloads each one using
# download_ncbi_genome.sh. Tracks which succeeded, which failed, resumes
# cleanly on re-run.
#
# Input format (config/ncbi_accessions.tsv):
#   species_id    accession           file_types
#   Cgar          GCF_024256425.1     fasta,gff,protein
#   Cmac          GCF_019173215.1     fasta,gff,protein
#   ...
#
# Output layout:
#   raw_genomes/
#     Cgar/
#       GCF_024256425.1_xxx_genomic.fna.gz
#       GCF_024256425.1_xxx_genomic.gff.gz
#       ...
#     Cmac/
#       ...
#     _download_log.tsv       # status table
#
# Usage:
#   sbatch STEP_00a_batch_download_genomes.sh [accessions_tsv] [output_dir]
#
# Re-run is idempotent: already-downloaded + md5-verified files are skipped.
# Failed downloads are retried on subsequent runs.
# ============================================================================
set -euo pipefail

ACC_TSV="${1:-config/ncbi_accessions.tsv}"
OUTDIR="${2:-raw_genomes}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNLOADER="${SCRIPT_DIR}/download_ncbi_genome.sh"

[[ -f "$ACC_TSV" ]] || { echo "ERROR: $ACC_TSV not found"; exit 1; }
[[ -f "$DOWNLOADER" ]] || { echo "ERROR: $DOWNLOADER not found"; exit 1; }

mkdir -p "$OUTDIR" logs
LOG="$OUTDIR/_download_log.tsv"

# Initialize log if it doesn't exist
if [[ ! -f "$LOG" ]]; then
    echo -e "species_id\taccession\tstatus\ttimestamp\tnotes" > "$LOG"
fi

N_OK=0
N_FAIL=0
N_SKIP=0

while IFS=$'\t' read -r sp_id acc file_types; do
    [[ "$sp_id" =~ ^# ]] && continue
    [[ -z "$sp_id" ]] && continue
    [[ "$sp_id" == "species_id" ]] && continue  # skip header

    sp_dir="$OUTDIR/$sp_id"
    mkdir -p "$sp_dir"

    # Check if we already have a successful log entry
    if grep -qP "^${sp_id}\t${acc}\tOK\t" "$LOG" 2>/dev/null; then
        # Verify files still exist
        if ls "$sp_dir"/${acc}_*_genomic.fna.gz &>/dev/null; then
            echo "[SKIP] $sp_id ($acc) — already downloaded"
            N_SKIP=$((N_SKIP+1))
            continue
        else
            echo "  Previous OK entry but files missing — re-downloading"
        fi
    fi

    echo ""
    echo "=============================================="
    echo "Downloading $sp_id ($acc)"
    echo "=============================================="

    ts=$(date -Iseconds)
    if bash "$DOWNLOADER" "$acc" "$sp_dir" "${file_types:-fasta,gff,protein}" 2>&1; then
        echo -e "${sp_id}\t${acc}\tOK\t${ts}\t" >> "$LOG"
        N_OK=$((N_OK+1))
    else
        echo -e "${sp_id}\t${acc}\tFAIL\t${ts}\tExit code $?" >> "$LOG"
        N_FAIL=$((N_FAIL+1))
    fi

done < "$ACC_TSV"

echo ""
echo "=============================================="
echo "Batch summary"
echo "=============================================="
echo "  OK:      $N_OK"
echo "  FAIL:    $N_FAIL"
echo "  SKIP:    $N_SKIP"
echo "  Log:     $LOG"
echo ""

if [[ $N_FAIL -gt 0 ]]; then
    echo "To retry failed downloads, just re-run this script."
    echo "Failed entries:"
    grep -P "\tFAIL\t" "$LOG" | tail -n "$N_FAIL"
    exit 2
fi
