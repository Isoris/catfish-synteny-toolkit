#!/usr/bin/env python3
"""
STEP_09b_score_dual_evidence.py

Consume mashmap multi-pi sweep outputs and miniprot GFF, compute per-tile
dual-evidence scores for each breakpoint window.

Tile scoring:
  KMER_ROBUSTNESS = highest pi at which this tile has an external homology
                    block. Scale: 0 (no homology even at pi=70) to 95 (homology
                    at pi=95). Higher = more conserved.
  PROTEIN_SUPPORT = 1 if this tile overlaps a miniprot mRNA feature, else 0.
  AGREEMENT       = "BOTH" if PROTEIN_SUPPORT=1 AND KMER_ROBUSTNESS>=80
                    "KMER_ONLY" if KMER_ROBUSTNESS>=80 and no protein
                    "PROTEIN_ONLY" if protein hit but weak kmer
                    "NEITHER" otherwise

Per-breakpoint confidence:
  Look at the flanks of the breakpoint (tiles left and right of midpoint):
    HIGH  : both flanks have majority BOTH or KMER_ONLY tiles
    MEDIUM: one flank clean, the other mixed
    LOW   : both flanks NEITHER-dominant (assembly gap? repeat desert?)
    NOVELTY: clean asymmetry — one flank orthologous elsewhere, other is
             lineage-specific (potentially real biological signal)
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ----------------------------------------------------------------------------
# Parse mashmap PAF (same format as wfmash)
# ----------------------------------------------------------------------------
def parse_paf_windows(path: Path) -> dict[str, list[tuple[int, int]]]:
    """
    Return {window_name: [(qstart, qend), ...]} of intervals on the query that
    have a homology hit to any other window.
    """
    hits: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if not path.exists():
        return hits
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname = f[0]
            tname = f[5]
            if qname == tname:
                continue  # self-hit
            qstart = int(f[2]); qend = int(f[3])
            # Only count hits to a different breakpoint window (different bp_id)
            qbp = qname.split("__")[0]
            tbp = tname.split("__")[0]
            if qbp == tbp:
                continue
            hits[qname].append((qstart, qend))
    return hits


# ----------------------------------------------------------------------------
# Parse miniprot GFF
# ----------------------------------------------------------------------------
def parse_miniprot_gff(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    """
    Return {window_name: [(start, end, protein_id), ...]} for mRNA features.
    """
    hits: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    if not path.exists():
        return hits
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] != "mRNA":
                continue
            win = f[0]
            start = int(f[3]) - 1  # GFF is 1-based closed; convert to 0-based half-open
            end = int(f[4])
            # Parse Target= from attributes
            attrs = f[8]
            protein_id = "unknown"
            for kv in attrs.split(";"):
                if kv.startswith("Target="):
                    protein_id = kv.split("=", 1)[1].split(" ")[0]
                    break
            hits[win].append((start, end, protein_id))
    return hits


# ----------------------------------------------------------------------------
# Tile evidence scoring
# ----------------------------------------------------------------------------
@dataclass
class TileScore:
    bp_id: str
    window: str
    tile_start: int
    tile_end: int
    kmer_robustness: int  # highest pi at which tile has hit (0 if none)
    protein_support: int  # 0 or 1
    n_proteins: int
    agreement: str


def interval_overlaps(a_start: int, a_end: int, intervals: list) -> bool:
    """Generic interval overlap check."""
    for iv in intervals:
        s, e = iv[0], iv[1]
        if s < a_end and e > a_start:
            return True
    return False


def count_interval_overlaps(a_start: int, a_end: int, intervals: list) -> int:
    n = 0
    for iv in intervals:
        s, e = iv[0], iv[1]
        if s < a_end and e > a_start:
            n += 1
    return n


def tile_windows(
    bp_bed: Path,
    kmer_hits_by_pi: dict[int, dict[str, list]],
    protein_hits: dict[str, list],
    tile_bp: int,
    pi_levels: list[int],
) -> list[TileScore]:
    scores = []

    with open(bp_bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            region = f[0]
            start = int(f[1])
            end = int(f[2])
            bp_id = f[3]

            # Reconstruct the window name used in mashmap/miniprot inputs
            win_name = f"{bp_id}__{region}_{start}_{end}"
            win_len = end - start

            for tile_start in range(0, win_len, tile_bp):
                tile_end = min(tile_start + tile_bp, win_len)

                # K-mer robustness: highest pi with a hit
                kmer_robustness = 0
                for pi in sorted(pi_levels, reverse=True):
                    intervals = kmer_hits_by_pi.get(pi, {}).get(win_name, [])
                    if interval_overlaps(tile_start, tile_end, intervals):
                        kmer_robustness = pi
                        break

                # Protein support
                prot_intervals = protein_hits.get(win_name, [])
                n_prot = count_interval_overlaps(tile_start, tile_end, prot_intervals)
                prot_support = 1 if n_prot > 0 else 0

                # Agreement category
                if prot_support and kmer_robustness >= 80:
                    agreement = "BOTH"
                elif kmer_robustness >= 80:
                    agreement = "KMER_ONLY"
                elif prot_support:
                    agreement = "PROTEIN_ONLY"
                else:
                    agreement = "NEITHER"

                scores.append(TileScore(
                    bp_id=bp_id, window=win_name,
                    tile_start=tile_start, tile_end=tile_end,
                    kmer_robustness=kmer_robustness,
                    protein_support=prot_support,
                    n_proteins=n_prot,
                    agreement=agreement,
                ))

    return scores


# ----------------------------------------------------------------------------
# Per-breakpoint confidence
# ----------------------------------------------------------------------------
def classify_breakpoint(tiles: list[TileScore]) -> dict:
    """
    Analyze tiles for one breakpoint. The breakpoint midpoint is at the center
    of the window. Compare left vs right flanks.
    """
    if not tiles:
        return {"confidence": "UNKNOWN", "note": "no tiles"}

    tiles = sorted(tiles, key=lambda t: t.tile_start)
    n = len(tiles)
    mid = n // 2
    left = tiles[:mid]
    right = tiles[mid:]

    def flank_profile(flank):
        counts = defaultdict(int)
        for t in flank:
            counts[t.agreement] += 1
        total = len(flank)
        if total == 0:
            return {"clean_frac": 0.0, "dominant": "NONE"}
        clean = counts["BOTH"] + counts["KMER_ONLY"]
        clean_frac = clean / total
        dominant = max(counts, key=counts.get)
        return {"clean_frac": clean_frac, "dominant": dominant,
                "n_both": counts["BOTH"], "n_kmer_only": counts["KMER_ONLY"],
                "n_protein_only": counts["PROTEIN_ONLY"], "n_neither": counts["NEITHER"],
                "total": total}

    L = flank_profile(left)
    R = flank_profile(right)

    # Classify
    if L["clean_frac"] >= 0.6 and R["clean_frac"] >= 0.6:
        confidence = "HIGH"
        note = "both flanks show robust dual-evidence homology"
    elif L["clean_frac"] >= 0.6 or R["clean_frac"] >= 0.6:
        if abs(L["clean_frac"] - R["clean_frac"]) >= 0.5:
            confidence = "NOVELTY"
            note = "asymmetric: one flank has homology, other lineage-specific"
        else:
            confidence = "MEDIUM"
            note = "one flank clean, other mixed"
    else:
        if L["dominant"] == "NEITHER" and R["dominant"] == "NEITHER":
            confidence = "LOW"
            note = "both flanks lack evidence — assembly gap, repeat desert, or extreme divergence"
        else:
            confidence = "MEDIUM"
            note = "partial evidence on both flanks"

    return {
        "confidence": confidence,
        "note": note,
        "left_clean_frac": L["clean_frac"],
        "right_clean_frac": R["clean_frac"],
        "left_dominant": L["dominant"],
        "right_dominant": R["dominant"],
        "n_tiles": n,
    }


# ----------------------------------------------------------------------------
# Visualization
# ----------------------------------------------------------------------------
def plot_evidence_tracks(
    tiles_by_bp: dict[str, list[TileScore]],
    out_pdf: Path,
    max_per_page: int = 12,
):
    bp_ids = sorted(tiles_by_bp.keys())
    n_bp = len(bp_ids)
    if n_bp == 0:
        return

    ncols = 2
    nrows = min(max_per_page // ncols, (n_bp + ncols - 1) // ncols)
    n_pages = (n_bp + max_per_page - 1) // max_per_page

    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(out_pdf) as pdf:
        for page_idx in range(n_pages):
            start_i = page_idx * max_per_page
            end_i = min(start_i + max_per_page, n_bp)
            page_bps = bp_ids[start_i:end_i]

            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(14, 2.5 * nrows),
                                     squeeze=False)
            for i, bp_id in enumerate(page_bps):
                ax = axes[i // ncols][i % ncols]
                tiles = sorted(tiles_by_bp[bp_id], key=lambda t: t.tile_start)

                x = [t.tile_start / 1000 for t in tiles]  # kb
                kmer = [t.kmer_robustness for t in tiles]
                prot = [t.n_proteins for t in tiles]

                ax.bar(x, kmer, width=4.8, color="steelblue",
                       alpha=0.6, label="k-mer pi max")
                ax2 = ax.twinx()
                ax2.bar([xi + 0.1 for xi in x], prot, width=4.8, color="orange",
                        alpha=0.7, label="# proteins")

                ax.axvline(max(x) / 2, color="red", linestyle="--", linewidth=1,
                           label="breakpoint")
                ax.set_title(bp_id, fontsize=9)
                ax.set_ylabel("k-mer pi", color="steelblue", fontsize=8)
                ax2.set_ylabel("proteins", color="orange", fontsize=8)
                ax.set_xlabel("position (kb)", fontsize=8)
                ax.set_ylim(0, 100)
                ax.tick_params(labelsize=7)
                ax2.tick_params(labelsize=7)

            # Blank unused
            for j in range(len(page_bps), nrows * ncols):
                axes[j // ncols][j % ncols].axis("off")

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bp-bed", required=True, type=Path)
    ap.add_argument("--mash-dir", required=True, type=Path)
    ap.add_argument("--miniprot-gff", required=True, type=Path)
    ap.add_argument("--pi-levels", nargs="+", type=int, required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--tile-bp", type=int, default=5000)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    print(f"[STEP_09b] Loading mashmap sweeps at pi={args.pi_levels}", file=sys.stderr)
    kmer_hits_by_pi: dict[int, dict[str, list]] = {}
    for pi in args.pi_levels:
        paf_path = args.mash_dir / f"windows_pi{pi}.paf"
        kmer_hits_by_pi[pi] = parse_paf_windows(paf_path)
        total = sum(len(v) for v in kmer_hits_by_pi[pi].values())
        print(f"  pi={pi}: {total} hits", file=sys.stderr)

    print(f"[STEP_09b] Loading miniprot GFF: {args.miniprot_gff}", file=sys.stderr)
    protein_hits = parse_miniprot_gff(args.miniprot_gff)
    total_prot = sum(len(v) for v in protein_hits.values())
    print(f"  {total_prot} protein alignments across {len(protein_hits)} windows", file=sys.stderr)

    print(f"[STEP_09b] Scoring tiles...", file=sys.stderr)
    scores = tile_windows(args.bp_bed, kmer_hits_by_pi, protein_hits,
                          args.tile_bp, args.pi_levels)

    # Group by bp_id
    tiles_by_bp: dict[str, list[TileScore]] = defaultdict(list)
    for s in scores:
        tiles_by_bp[s.bp_id].append(s)

    # Write tile TSV
    tile_tsv = args.out / "tile_evidence.tsv"
    with open(tile_tsv, "w") as fh:
        fh.write("bp_id\twindow\ttile_start\ttile_end\tkmer_robustness\t"
                 "protein_support\tn_proteins\tagreement\n")
        for s in scores:
            fh.write(f"{s.bp_id}\t{s.window}\t{s.tile_start}\t{s.tile_end}\t"
                     f"{s.kmer_robustness}\t{s.protein_support}\t{s.n_proteins}\t{s.agreement}\n")
    print(f"[STEP_09b] Wrote {len(scores)} tile rows -> {tile_tsv}", file=sys.stderr)

    # Per-bp confidence
    conf_tsv = args.out / "breakpoint_confidence.tsv"
    with open(conf_tsv, "w") as fh:
        fh.write("bp_id\tconfidence\tnote\tleft_clean_frac\tright_clean_frac\t"
                 "left_dominant\tright_dominant\tn_tiles\n")
        conf_counts = defaultdict(int)
        for bp_id, tiles in sorted(tiles_by_bp.items()):
            cls = classify_breakpoint(tiles)
            conf_counts[cls["confidence"]] += 1
            fh.write(f"{bp_id}\t{cls['confidence']}\t{cls['note']}\t"
                     f"{cls['left_clean_frac']:.2f}\t{cls['right_clean_frac']:.2f}\t"
                     f"{cls['left_dominant']}\t{cls['right_dominant']}\t"
                     f"{cls['n_tiles']}\n")
    print(f"[STEP_09b] Per-breakpoint confidence -> {conf_tsv}", file=sys.stderr)
    print(f"\n[STEP_09b] Confidence summary:", file=sys.stderr)
    for conf, n in sorted(conf_counts.items()):
        print(f"  {conf:10s} {n:4d}", file=sys.stderr)

    # Visualize
    pdf_out = args.out / "evidence_tracks.pdf"
    print(f"[STEP_09b] Plotting evidence tracks -> {pdf_out}", file=sys.stderr)
    plot_evidence_tracks(tiles_by_bp, pdf_out)

    print(f"\n[STEP_09b] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
