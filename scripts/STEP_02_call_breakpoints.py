#!/usr/bin/env python3
"""
STEP_02_call_breakpoints.py

Parse wfmash scaffold PAF to call structural breakpoints between species pairs.
Detects three event classes:

    FUSION   : in query species, two adjacent blocks map to different target chroms
               (or vice versa — fission when viewed from the other direction)
    FISSION  : symmetric of fusion, detected by swapping query/target roles
    INVERSION: adjacent collinear blocks on the same target chrom with strand flip

Output BED-like TSV with columns:
    query_species  query_chrom  q_start  q_end
    target_species target_chrom t_start  t_end
    event_type    left_block_id  right_block_id  gap_bp

Usage:
    python3 STEP_02_call_breakpoints.py \
        --paf results/01_wfmash/catfish_core.scaffolds.paf \
        --out results/02_breakpoints/
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


# ----------------------------------------------------------------------------
# PAF parsing
# ----------------------------------------------------------------------------
@dataclass
class PafRecord:
    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    matches: int
    aln_len: int
    mapq: int

    @property
    def qspecies(self) -> str:
        # PanSN: species#hap#chrom
        return self.qname.split("#", 1)[0]

    @property
    def tspecies(self) -> str:
        return self.tname.split("#", 1)[0]

    @property
    def qchrom(self) -> str:
        parts = self.qname.split("#", 2)
        return parts[2] if len(parts) == 3 else self.qname

    @property
    def tchrom(self) -> str:
        parts = self.tname.split("#", 2)
        return parts[2] if len(parts) == 3 else self.tname


def parse_paf(path: Path) -> list[PafRecord]:
    records = []
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            records.append(PafRecord(
                qname=f[0], qlen=int(f[1]), qstart=int(f[2]), qend=int(f[3]),
                strand=f[4],
                tname=f[5], tlen=int(f[6]), tstart=int(f[7]), tend=int(f[8]),
                matches=int(f[9]), aln_len=int(f[10]), mapq=int(f[11]),
            ))
    return records


# ----------------------------------------------------------------------------
# Breakpoint detection
# ----------------------------------------------------------------------------
def call_breakpoints_for_pair(
    blocks: list[PafRecord],
    min_block_bp: int = 200_000,
    max_gap_bp: int = 2_000_000,
) -> list[dict]:
    """
    Given all PAF blocks between ONE query species and ONE target species,
    walk each query chromosome and detect breakpoints between consecutive blocks.
    """
    # Filter by minimum block size
    blocks = [b for b in blocks if (b.qend - b.qstart) >= min_block_bp]

    # Group by query chromosome
    by_qchrom: dict[str, list[PafRecord]] = defaultdict(list)
    for b in blocks:
        by_qchrom[b.qchrom].append(b)

    events = []
    for qchrom, blist in by_qchrom.items():
        blist.sort(key=lambda b: b.qstart)
        for b1, b2 in zip(blist[:-1], blist[1:]):
            gap = b2.qstart - b1.qend
            if gap > max_gap_bp:
                # Too far apart to call a breakpoint reliably
                continue

            qspecies = b1.qspecies
            tspecies = b1.tspecies

            # CASE 1: FUSION — adjacent blocks on same query chrom map to different target chroms
            if b1.tchrom != b2.tchrom:
                events.append({
                    "query_species": qspecies, "query_chrom": qchrom,
                    "q_start": b1.qend, "q_end": b2.qstart,
                    "target_species": tspecies,
                    "left_target_chrom": b1.tchrom,
                    "left_target_pos": b1.tend if b1.strand == "+" else b1.tstart,
                    "right_target_chrom": b2.tchrom,
                    "right_target_pos": b2.tstart if b2.strand == "+" else b2.tend,
                    "event_type": "FUSION_in_query",
                    "gap_bp": gap,
                    "left_strand": b1.strand,
                    "right_strand": b2.strand,
                })

            # CASE 2: INVERSION — same target chrom, strand flip between adjacent blocks
            elif b1.strand != b2.strand:
                events.append({
                    "query_species": qspecies, "query_chrom": qchrom,
                    "q_start": b1.qend, "q_end": b2.qstart,
                    "target_species": tspecies,
                    "left_target_chrom": b1.tchrom,
                    "left_target_pos": b1.tend if b1.strand == "+" else b1.tstart,
                    "right_target_chrom": b2.tchrom,
                    "right_target_pos": b2.tstart if b2.strand == "+" else b2.tend,
                    "event_type": "INVERSION_breakpoint",
                    "gap_bp": gap,
                    "left_strand": b1.strand,
                    "right_strand": b2.strand,
                })

            # CASE 3: same chrom, same strand, but non-collinear target positions
            # (could be a small inversion we missed, or a translocation)
            elif b1.strand == "+" and b2.tstart < b1.tend:
                events.append({
                    "query_species": qspecies, "query_chrom": qchrom,
                    "q_start": b1.qend, "q_end": b2.qstart,
                    "target_species": tspecies,
                    "left_target_chrom": b1.tchrom,
                    "left_target_pos": b1.tend,
                    "right_target_chrom": b2.tchrom,
                    "right_target_pos": b2.tstart,
                    "event_type": "NONCOLLINEAR",
                    "gap_bp": gap,
                    "left_strand": b1.strand,
                    "right_strand": b2.strand,
                })

    return events


def call_all_breakpoints(
    records: list[PafRecord],
    min_block_bp: int = 200_000,
    max_gap_bp: int = 2_000_000,
) -> list[dict]:
    """Detect breakpoints for every (query_species, target_species) pair."""
    # Group by (qspecies, tspecies)
    by_pair: dict[tuple[str, str], list[PafRecord]] = defaultdict(list)
    for r in records:
        if r.qspecies == r.tspecies:
            continue  # skip intra-species (shouldn't happen with -Y '#' but belt-and-suspenders)
        by_pair[(r.qspecies, r.tspecies)].append(r)

    all_events = []
    for (qsp, tsp), blist in by_pair.items():
        print(f"  Pair {qsp} -> {tsp}: {len(blist)} scaffold blocks", file=sys.stderr)
        events = call_breakpoints_for_pair(blist, min_block_bp, max_gap_bp)
        all_events.extend(events)

    return all_events


# ----------------------------------------------------------------------------
# Reciprocal fusion/fission classification
# ----------------------------------------------------------------------------
def classify_fusion_fission(events: list[dict]) -> list[dict]:
    """
    A 'FUSION_in_query' event in species A -> species B direction is a
    'FISSION_in_query' in the reverse direction. We can reconcile by looking
    for reciprocal events.

    Simpler heuristic here: just annotate the reciprocal expectation.
    True polarization requires an outgroup — that's a separate step.
    """
    for e in events:
        if e["event_type"] == "FUSION_in_query":
            e["note"] = (
                f"Chrom {e['query_chrom']} in {e['query_species']} spans "
                f"two distinct chroms ({e['left_target_chrom']}, {e['right_target_chrom']}) "
                f"in {e['target_species']}. Either fusion in {e['query_species']} "
                f"or fission in {e['target_species']}. Polarize with outgroup."
            )
        elif e["event_type"] == "INVERSION_breakpoint":
            e["note"] = "Adjacent collinear blocks with strand flip"
        else:
            e["note"] = ""
    return events


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True, type=Path, help="wfmash scaffold PAF file")
    ap.add_argument("--out", required=True, type=Path, help="Output directory")
    ap.add_argument("--min-block-bp", type=int, default=200_000,
                    help="Minimum block size to consider (default 200 kb)")
    ap.add_argument("--max-gap-bp", type=int, default=2_000_000,
                    help="Max gap between adjacent blocks to call a breakpoint (default 2 Mb)")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    print(f"[STEP_02] Parsing {args.paf}", file=sys.stderr)
    records = parse_paf(args.paf)
    print(f"[STEP_02] Loaded {len(records)} scaffold records", file=sys.stderr)

    print(f"[STEP_02] Calling breakpoints...", file=sys.stderr)
    events = call_all_breakpoints(records, args.min_block_bp, args.max_gap_bp)
    events = classify_fusion_fission(events)

    # Write main TSV
    out_tsv = args.out / "breakpoints_all.tsv"
    if events:
        with open(out_tsv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(events[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(events)
    print(f"[STEP_02] Wrote {len(events)} events -> {out_tsv}", file=sys.stderr)

    # Split by event type
    by_type: dict[str, list] = defaultdict(list)
    for e in events:
        by_type[e["event_type"]].append(e)

    print("\n[STEP_02] Event summary:", file=sys.stderr)
    for etype, elist in sorted(by_type.items()):
        print(f"  {etype:30s} {len(elist):6d}", file=sys.stderr)
        out_sub = args.out / f"breakpoints_{etype}.tsv"
        with open(out_sub, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(elist[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(elist)

    # Per-species-pair summary
    summary_path = args.out / "breakpoint_summary.tsv"
    pair_counts: dict[tuple, dict] = defaultdict(lambda: defaultdict(int))
    for e in events:
        key = (e["query_species"], e["target_species"])
        pair_counts[key][e["event_type"]] += 1
        pair_counts[key]["TOTAL"] += 1

    with open(summary_path, "w") as fh:
        fh.write("query_species\ttarget_species\tfusion_events\tinversion_events\tnoncollinear\ttotal\n")
        for (q, t), counts in sorted(pair_counts.items()):
            fh.write(f"{q}\t{t}\t"
                     f"{counts['FUSION_in_query']}\t"
                     f"{counts['INVERSION_breakpoint']}\t"
                     f"{counts['NONCOLLINEAR']}\t"
                     f"{counts['TOTAL']}\n")
    print(f"[STEP_02] Per-pair summary -> {summary_path}", file=sys.stderr)

    # Also emit BED files for downstream local refinement with minimap2
    bed_path = args.out / "breakpoint_windows.bed"
    with open(bed_path, "w") as fh:
        for e in events:
            # ±500 kb window around breakpoint for local re-alignment
            pad = 500_000
            center = (e["q_start"] + e["q_end"]) // 2
            fh.write(f"{e['query_species']}#1#{e['query_chrom']}\t"
                     f"{max(0, center - pad)}\t{center + pad}\t"
                     f"{e['event_type']}_{e['target_species']}_{e['left_target_chrom']}-{e['right_target_chrom']}\n")
    print(f"[STEP_02] Windows for minimap2 refinement -> {bed_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
