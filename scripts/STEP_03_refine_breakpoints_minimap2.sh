#!/usr/bin/env bash
#SBATCH --job-name=bp_refine
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/bp_refine_%j.out
#SBATCH --error=logs/bp_refine_%j.err
# ============================================================================
# STEP_03_refine_breakpoints_minimap2.sh
#   For each candidate breakpoint window from STEP_02, extract ±500 kb from
#   both query and target genomes and realign with minimap2 asm20 for
#   base-pair resolution of the breakpoint.
#
# Usage:
#   sbatch STEP_03_refine_breakpoints_minimap2.sh <tier>
# ============================================================================
set -euo pipefail

TIER="${1:-core}"
module load Miniconda3 || true
source activate assembly

BED="results/02_breakpoints/breakpoint_windows.bed"
CONCAT_FA="results/00_panSN_inputs/catfish_${TIER}.fa.gz"
OUTDIR="results/03_refined_breakpoints"
mkdir -p "$OUTDIR"

if [[ ! -f "$BED" ]]; then
    echo "ERROR: $BED not found. Run STEP_02 first."
    exit 1
fi

echo "[STEP_03] Refining $(wc -l < "$BED") breakpoint windows"

# For each breakpoint window, extract ±500 kb from query side AND the two
# candidate target regions, then realign.
python3 - "$BED" "$CONCAT_FA" "$OUTDIR" "${SLURM_CPUS_PER_TASK:-32}" <<'PYEOF'
import sys
import subprocess
from pathlib import Path

bed_path    = Path(sys.argv[1])
fa_path     = Path(sys.argv[2])
outdir      = Path(sys.argv[3])
threads     = sys.argv[4]

with open(bed_path) as fh:
    for i, line in enumerate(fh):
        f = line.rstrip("\n").split("\t")
        q_region = f[0]           # e.g. Cgar#1#LG01
        q_start  = int(f[1])
        q_end    = int(f[2])
        label    = f[3]

        region_str = f"{q_region}:{q_start}-{q_end}"
        q_fa = outdir / f"bp_{i:04d}_query.fa"
        t_fa = outdir / f"bp_{i:04d}_target.fa"
        out_paf = outdir / f"bp_{i:04d}_refined.paf"

        # Extract query window
        subprocess.run(
            ["samtools", "faidx", str(fa_path), region_str, "-o", str(q_fa)],
            check=True,
        )

        # Parse target info from label: EVENT_TSPECIES_LEFTCHROM-RIGHTCHROM
        # label format: FUSION_in_query_Cmac_LG05-LG12
        parts = label.rsplit("_", 2)
        if len(parts) < 3:
            continue
        tspecies = parts[-2]
        chroms = parts[-1].split("-")

        # Extract full target chroms (simpler than trying to guess positions)
        target_regions = [f"{tspecies}#1#{c}" for c in chroms if c]
        with open(t_fa, "w") as tfh:
            for tr in target_regions:
                result = subprocess.run(
                    ["samtools", "faidx", str(fa_path), tr],
                    capture_output=True, text=True,
                )
                if result.returncode == 0:
                    tfh.write(result.stdout)

        # Realign: asm20 preset for ~90%+ identity
        with open(out_paf, "w") as ofh:
            subprocess.run(
                ["minimap2", "-cx", "asm20", "-t", threads,
                 str(t_fa), str(q_fa)],
                stdout=ofh, check=True,
            )

        print(f"  [{i:04d}] {label}: refined -> {out_paf}")
PYEOF

# Merge all refined PAFs
cat "$OUTDIR"/bp_*_refined.paf > "$OUTDIR/all_refined.paf"
echo "[STEP_03] Merged PAF: $OUTDIR/all_refined.paf ($(wc -l < "$OUTDIR/all_refined.paf") records)"
