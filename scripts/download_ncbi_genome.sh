#!/usr/bin/env bash
# ============================================================================
# download_ncbi_genome.sh
#
# Robust downloader for a single NCBI assembly accession.
# Uses deterministic FTP URL construction + wget with aggressive retries
# + md5 verification.
#
# Does NOT use the `datasets` CLI (unreliable for batch).
# Does NOT use rsync (deprecated by NCBI June 2026).
#
# Usage:
#   bash download_ncbi_genome.sh <accession> <output_dir> [file_types]
#
#   accession:  GCA_XXXXXXXXX.Y or GCF_XXXXXXXXX.Y
#   output_dir: where files will land (will be created)
#   file_types: comma-separated subset of:
#                 fasta    (*_genomic.fna.gz)         DEFAULT
#                 gff      (*_genomic.gff.gz)
#                 protein  (*_protein.faa.gz)
#                 cds      (*_cds_from_genomic.fna.gz)
#                 all      (all of the above)
#               Default: "fasta,gff,protein"
#
# Exit codes:
#   0 = success, all requested files downloaded and md5 verified
#   1 = assembly directory not found on FTP (bad accession? unreleased?)
#   2 = download failed after retries
#   3 = md5 mismatch after download
# ============================================================================
set -euo pipefail

ACC="${1:?Usage: download_ncbi_genome.sh <accession> <output_dir> [file_types]}"
OUTDIR="${2:?Usage: download_ncbi_genome.sh <accession> <output_dir> [file_types]}"
FILE_TYPES="${3:-fasta,gff,protein}"

mkdir -p "$OUTDIR"

# ---- Expand "all" ----
if [[ "$FILE_TYPES" == "all" ]]; then
    FILE_TYPES="fasta,gff,protein,cds"
fi

# ---- Parse accession ----
# Example: GCA_003123456.1 -> prefix=GCA, digits=003123456
if [[ ! "$ACC" =~ ^(GCA|GCF)_([0-9]+)\.[0-9]+$ ]]; then
    echo "ERROR: Accession must look like GCA_123456789.1 or GCF_123456789.1"
    exit 1
fi
PREFIX="${BASH_REMATCH[1]}"
DIGITS="${BASH_REMATCH[2]}"

# ---- Construct FTP path: split digits into 3/3/3 ----
# 003123456 -> 003/123/456
PART1="${DIGITS:0:3}"
PART2="${DIGITS:3:3}"
PART3="${DIGITS:6:3}"

BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/${PREFIX}/${PART1}/${PART2}/${PART3}"

echo "[download] accession: $ACC"
echo "[download] base URL:  $BASE_URL"

# ---- Get directory listing to find the actual assembly dir (has assembly name suffix) ----
LISTING=$(mktemp)
trap "rm -f $LISTING" EXIT

# retry the listing fetch up to 5 times
ok=0
for attempt in 1 2 3 4 5; do
    if curl -sS --retry 3 --retry-delay 10 --max-time 60 \
           "${BASE_URL}/" -o "$LISTING" 2>/dev/null; then
        ok=1; break
    fi
    echo "  Listing attempt $attempt failed, retrying in 15s..."
    sleep 15
done

if [[ $ok -ne 1 ]] || [[ ! -s "$LISTING" ]]; then
    echo "ERROR: Could not list ${BASE_URL}/ — assembly may not exist or network blocked"
    exit 1
fi

# Find the assembly subdirectory (pattern: GCA_003123456.1_<assembly_name>)
ASM_DIR=$(grep -oE "${ACC}_[^/\"<>]+" "$LISTING" | head -1)
if [[ -z "$ASM_DIR" ]]; then
    echo "ERROR: Could not locate assembly subdirectory for $ACC in $BASE_URL"
    echo "Listing snippet:"
    head -20 "$LISTING"
    exit 1
fi

FULL_URL="${BASE_URL}/${ASM_DIR}"
echo "[download] assembly dir: $ASM_DIR"

# ---- Download md5 checksums first ----
MD5_FILE="$OUTDIR/${ACC}.md5checksums.txt"
if ! wget --tries=5 --waitretry=20 --timeout=60 --continue --quiet \
        -O "$MD5_FILE" "${FULL_URL}/md5checksums.txt"; then
    echo "  WARNING: could not fetch md5checksums.txt (continuing without verification)"
    : > "$MD5_FILE"
fi

# ---- Map file_types to NCBI suffixes ----
declare -A TYPE_SUFFIX=(
    [fasta]="_genomic.fna.gz"
    [gff]="_genomic.gff.gz"
    [protein]="_protein.faa.gz"
    [cds]="_cds_from_genomic.fna.gz"
)

# ---- Download each requested file with retry + resume ----
status=0
IFS=',' read -ra REQ <<< "$FILE_TYPES"
for t in "${REQ[@]}"; do
    t="${t// /}"
    suffix="${TYPE_SUFFIX[$t]:-}"
    if [[ -z "$suffix" ]]; then
        echo "  WARNING: unknown file type '$t', skipping"
        continue
    fi

    fname="${ASM_DIR}${suffix}"
    url="${FULL_URL}/${fname}"
    out="$OUTDIR/${fname}"

    if [[ -s "$out" ]]; then
        # Already have it — verify if possible
        if [[ -s "$MD5_FILE" ]]; then
            expected=$(grep "  \./${fname}" "$MD5_FILE" | awk '{print $1}' || true)
            actual=$(md5sum "$out" | awk '{print $1}')
            if [[ -n "$expected" && "$expected" == "$actual" ]]; then
                echo "  [cached] $t: $out (md5 ok)"
                continue
            else
                echo "  [stale]  $t: md5 mismatch, re-downloading"
                rm -f "$out"
            fi
        else
            echo "  [cached] $t: $out (no md5 to verify)"
            continue
        fi
    fi

    echo "  [fetch]  $t: $url"

    # wget flags:
    #   --tries=10         retry up to 10 times
    #   --waitretry=30     wait 30s between retries (exponential-ish)
    #   --timeout=120      2-minute timeout per connection
    #   --continue         resume partial downloads
    #   --read-timeout=30  kill stalled reads
    if ! wget --tries=10 --waitretry=30 --timeout=120 \
              --read-timeout=30 --continue --quiet \
              -O "$out" "$url"; then
        echo "    DOWNLOAD FAILED: $t"
        status=2
        continue
    fi

    # Verify md5 if available
    if [[ -s "$MD5_FILE" ]]; then
        expected=$(grep "  \./${fname}" "$MD5_FILE" | awk '{print $1}' || true)
        if [[ -n "$expected" ]]; then
            actual=$(md5sum "$out" | awk '{print $1}')
            if [[ "$expected" != "$actual" ]]; then
                echo "    MD5 MISMATCH: expected $expected, got $actual"
                rm -f "$out"
                status=3
                continue
            fi
            echo "    [ok]     md5 verified"
        fi
    fi
done

if [[ $status -eq 0 ]]; then
    echo "[download] $ACC: all requested files downloaded"
else
    echo "[download] $ACC: completed with errors (status=$status)"
fi
exit $status
