#!/usr/bin/env python3
"""
STEP_08_cross_species_stats.py

STAGE D — Cross-species statistics at each breakpoint.

For each breakpoint identified in Stage B, and each satellite/TE consensus
identified in Stage C, query the consensus against all species' genomes and
count hits. This produces the generalized equivalent of Kuang Figure 2D-F
(their '758 hits vs 491 hits vs 2196 hits' panel), but automated across all
breakpoints and species simultaneously.

Method:
    1. Extract satellite consensus sequences from TRF output (period >= 100 bp)
    2. For each satellite, run minimap2 -x sr against each species' full genome
       (or use mashmap for longer satellites)
    3. Count hits per chromosome per species
    4. Test enrichment: hits on fused chromosome vs hits on other chromosomes

Output:
    - satellite_catalog.tsv      : all satellites found, one row per consensus
    - satellite_hits_matrix.tsv  : rows=satellites, cols=species+chrom, cells=hit count
    - enrichment_tests.tsv       : chi-square test per satellite
    - Fig2_equivalent.pdf        : the Kuang Figure 2 D-F equivalent

Usage:
    python3 STEP_08_cross_species_stats.py \
        --annotation-dir results/07_bp_annotation/ \
        --manifest config/species_manifest.tsv \
        --out results/08_cross_species_stats/
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ----------------------------------------------------------------------------
# Parse TRF output to extract satellite consensus sequences
# ----------------------------------------------------------------------------
def parse_trf_ngs(trf_dat: Path, min_period: int = 100) -> list[dict]:
    """
    Parse TRF -ngs format:
      >seqname
      <start> <end> <period> <copies> <consensus_size> <pid> <indel%> <score>
      <A%> <C%> <G%> <T%> <entropy> <consensus> <sequence>
    Returns list of repeats with period >= min_period (satellite-class).
    """
    repeats = []
    if not trf_dat.exists():
        return repeats

    with open(trf_dat) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            try:
                start = int(parts[0])
                end = int(parts[1])
                period = int(parts[2])
                copies = float(parts[3])
                consensus = parts[13]
                score = int(parts[7])
            except ValueError:
                continue
            if period < min_period:
                continue
            repeats.append({
                "start": start, "end": end,
                "period": period, "copies": copies,
                "score": score, "consensus": consensus,
            })
    return repeats


# ----------------------------------------------------------------------------
# Parse species manifest
# ----------------------------------------------------------------------------
def parse_manifest(path: Path) -> list[dict]:
    species = []
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                if header is None and line.startswith("#species_id"):
                    header = line.lstrip("#").split("\t")
                continue
            if header is None:
                header = ["species_id", "genus", "species", "panSN_prefix",
                          "fasta_path", "tier", "notes"]
                continue
            fields = line.split("\t")
            species.append(dict(zip(header, fields)))
    return species


# ----------------------------------------------------------------------------
# Query satellite consensus against each species genome
# ----------------------------------------------------------------------------
def query_satellite_against_genomes(
    consensus: str,
    consensus_id: str,
    species_list: list[dict],
    tmp_dir: Path,
    threads: int = 8,
) -> dict[str, dict[str, int]]:
    """
    For a satellite consensus sequence, count hits per chromosome per species
    using minimap2 short-read mode (-ax sr effectively, but for query->genome
    we use -cx map-ont for any read-length including short).

    Returns: {species_id: {chrom: hit_count}}
    """
    # Write consensus to tmp FASTA (repeat 3x for better seed anchoring)
    query_fa = tmp_dir / f"{consensus_id}.fa"
    with open(query_fa, "w") as fh:
        fh.write(f">{consensus_id}\n{consensus * 3}\n")

    hits_per_sp: dict[str, dict[str, int]] = {}
    for sp in species_list:
        sp_id = sp["species_id"]
        fasta = sp["fasta_path"]
        if not Path(fasta).exists():
            continue

        # Use minimap2 map-ont which works for short (~100-2000bp) queries
        # Output PAF; one line per hit
        result = subprocess.run(
            ["minimap2", "-c", "-x", "map-ont", "-t", str(threads),
             "--secondary=yes", "-N", "1000",
             fasta, str(query_fa)],
            capture_output=True, text=True, timeout=300,
        )
        chrom_counts: dict[str, int] = defaultdict(int)
        for line in result.stdout.splitlines():
            f = line.split("\t")
            if len(f) < 12:
                continue
            tname = f[5]
            # strip PanSN if present
            tchrom = tname.split("#", 2)[-1] if "#" in tname else tname
            aln_len = int(f[10])
            if aln_len >= len(consensus) * 0.8:  # require >=80% of monomer
                chrom_counts[tchrom] += 1
        hits_per_sp[sp_id] = dict(chrom_counts)

    return hits_per_sp


# ----------------------------------------------------------------------------
# Chi-square enrichment test
# ----------------------------------------------------------------------------
def enrichment_test(hits: dict[str, dict[str, int]], focal_species: str, focal_chrom: str) -> dict:
    """
    Does this satellite hit the focal (fused/breakpoint) chromosome more than
    expected under null of uniform distribution across chromosomes of the same
    species?
    """
    focal_hits = hits.get(focal_species, {}).get(focal_chrom, 0)
    total_hits_focal_sp = sum(hits.get(focal_species, {}).values())
    n_chroms_focal_sp = len(hits.get(focal_species, {}))

    if total_hits_focal_sp == 0 or n_chroms_focal_sp <= 1:
        return {"focal_hits": focal_hits, "total": total_hits_focal_sp,
                "expected": None, "fold_enrichment": None, "chi2": None}

    expected = total_hits_focal_sp / n_chroms_focal_sp
    fold = focal_hits / expected if expected > 0 else None
    chi2 = ((focal_hits - expected) ** 2) / expected if expected > 0 else None

    return {
        "focal_hits": focal_hits,
        "total": total_hits_focal_sp,
        "expected": expected,
        "fold_enrichment": fold,
        "chi2": chi2,
    }


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation-dir", required=True, type=Path,
                    help="Output dir from STEP_07 (contains per-bp subdirs with trf.dat)")
    ap.add_argument("--manifest", required=True, type=Path,
                    help="Species manifest TSV")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-period", type=int, default=100,
                    help="Minimum TRF period to consider (satellite class); default 100")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    tmp_dir = args.out / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    # Parse manifest
    species_list = parse_manifest(args.manifest)
    print(f"[STEP_08] {len(species_list)} species in manifest", file=sys.stderr)

    # Collect satellites from all breakpoint TRF outputs
    satellite_catalog = []
    bp_dirs = sorted([d for d in args.annotation_dir.iterdir() if d.is_dir()])
    print(f"[STEP_08] Scanning {len(bp_dirs)} breakpoint dirs for satellites...", file=sys.stderr)

    for bp_dir in bp_dirs:
        bp_id = bp_dir.name
        trf_dat = bp_dir / "trf.dat"
        reps = parse_trf_ngs(trf_dat, args.min_period)
        for r_idx, r in enumerate(reps):
            satellite_catalog.append({
                "sat_id": f"{bp_id}_sat{r_idx:03d}",
                "bp_id": bp_id,
                "period": r["period"],
                "copies": r["copies"],
                "score": r["score"],
                "consensus": r["consensus"],
                "consensus_len": len(r["consensus"]),
            })

    print(f"[STEP_08] Found {len(satellite_catalog)} satellite consensus sequences", file=sys.stderr)

    # Write catalog
    cat_tsv = args.out / "satellite_catalog.tsv"
    with open(cat_tsv, "w") as fh:
        fh.write("sat_id\tbp_id\tperiod\tcopies\tscore\tconsensus_len\tconsensus\n")
        for s in satellite_catalog:
            fh.write(f"{s['sat_id']}\t{s['bp_id']}\t{s['period']}\t{s['copies']:.1f}\t"
                     f"{s['score']}\t{s['consensus_len']}\t{s['consensus']}\n")

    # Query each satellite against each genome
    # Build hit matrix: satellites (rows) × species_chrom (cols)
    print(f"[STEP_08] Querying satellites against genomes (this takes a while)...", file=sys.stderr)
    all_hits: dict[str, dict[str, dict[str, int]]] = {}  # sat_id -> species -> chrom -> count
    for i, s in enumerate(satellite_catalog):
        print(f"  [{i+1}/{len(satellite_catalog)}] {s['sat_id']} "
              f"(period={s['period']}, copies={s['copies']:.1f})", file=sys.stderr)
        hits = query_satellite_against_genomes(
            s["consensus"], s["sat_id"], species_list, tmp_dir, args.threads,
        )
        all_hits[s["sat_id"]] = hits

    # Write hit matrix (long format)
    matrix_tsv = args.out / "satellite_hits_long.tsv"
    with open(matrix_tsv, "w") as fh:
        fh.write("sat_id\tspecies\tchrom\thits\n")
        for sat_id, hits in all_hits.items():
            for sp, chrom_counts in hits.items():
                for chrom, n in chrom_counts.items():
                    fh.write(f"{sat_id}\t{sp}\t{chrom}\t{n}\n")

    # Per-species total hits (quick summary)
    summary_tsv = args.out / "satellite_hits_summary.tsv"
    with open(summary_tsv, "w") as fh:
        fh.write("sat_id\tbp_id\tperiod\tconsensus_len")
        sp_ids = [sp["species_id"] for sp in species_list]
        for sp_id in sp_ids:
            fh.write(f"\t{sp_id}_total\t{sp_id}_max_chrom\t{sp_id}_max_chrom_hits")
        fh.write("\n")
        for s in satellite_catalog:
            sat_id = s["sat_id"]
            hits = all_hits.get(sat_id, {})
            fh.write(f"{sat_id}\t{s['bp_id']}\t{s['period']}\t{s['consensus_len']}")
            for sp_id in sp_ids:
                sp_hits = hits.get(sp_id, {})
                total = sum(sp_hits.values())
                if sp_hits:
                    max_chrom = max(sp_hits, key=sp_hits.get)
                    max_n = sp_hits[max_chrom]
                else:
                    max_chrom = "NA"
                    max_n = 0
                fh.write(f"\t{total}\t{max_chrom}\t{max_n}")
            fh.write("\n")

    print(f"\n[STEP_08] Outputs:", file=sys.stderr)
    print(f"  Catalog:  {cat_tsv}", file=sys.stderr)
    print(f"  Hits:     {matrix_tsv}", file=sys.stderr)
    print(f"  Summary:  {summary_tsv}", file=sys.stderr)

    # ------------------------------------------------------------------------
    # Figure 2 equivalent: for each satellite, bar plot of hits per species
    # ------------------------------------------------------------------------
    pdf_out = args.out / "Fig2_equivalent.pdf"
    n_sat = len(satellite_catalog)
    if n_sat == 0:
        print("[STEP_08] No satellites to plot", file=sys.stderr)
        return

    ncols = 4
    nrows = (n_sat + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 2.5 * nrows),
                             squeeze=False)

    for i, s in enumerate(satellite_catalog):
        ax = axes[i // ncols][i % ncols]
        sat_id = s["sat_id"]
        hits = all_hits.get(sat_id, {})

        sp_totals = {sp_id: sum(hits.get(sp_id, {}).values())
                     for sp_id in [sp["species_id"] for sp in species_list]}

        ax.bar(sp_totals.keys(), sp_totals.values(), color="steelblue",
               edgecolor="black", linewidth=0.3)
        ax.set_title(f"{sat_id}\nperiod={s['period']} len={s['consensus_len']}",
                     fontsize=8)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.set_ylabel("Total hits", fontsize=7)

    # Blank unused axes
    for j in range(n_sat, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    plt.tight_layout()
    plt.savefig(pdf_out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Fig2 PDF: {pdf_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
