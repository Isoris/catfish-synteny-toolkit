#!/usr/bin/env python3
"""
STEP_09c_flank_coherence.py

COLLECT-THEN-REFINE flank coherence scoring for breakpoints.

Adapts hamburger's (djw533) clustering logic to breakpoint validation:

  COLLECT (lenient): For each species, miniprot with --outs=0.5 gives all
                     protein family hits across the whole genome.
                     (Done in STEP_09c1, feeds into this step.)

  EXTRACT:           For each breakpoint flank, pull the ordered list of
                     protein families present in that flank.

  REFINE (strict):   Hamburger-style check:
                     1. Each flank must contain >= min_families protein
                        families that appear in the same ORDER in at least
                        one other species' genome.
                     2. Those families must be clustered (no more than
                        max_gap_genes intervening non-matching genes).

Output: per-flank coherence score + per-breakpoint confidence.

This is a 3rd evidence layer on top of STEP_09b's k-mer + per-protein hits.
The power comes from requiring ORDER preservation across species, which
random false positives cannot produce.

Usage:
    python3 STEP_09c_flank_coherence.py \
        --bp-bed results/05_synteny_graph/breakpoints.bed \
        --wg-miniprot-dir results/09c_wg_miniprot/ \
        --manifest config/species_manifest.tsv \
        --out results/09c_flank_coherence/ \
        --flank-kb 100 \
        --min-families 3 \
        --max-gap-genes 10
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


# ============================================================================
# Data structures
# ============================================================================
@dataclass
class GeneHit:
    species: str
    chrom: str            # without PanSN prefix
    start: int
    end: int
    strand: str
    family: str           # the Target= protein family ID from the reference
    score: float


@dataclass
class Flank:
    bp_id: str
    side: str             # 'left' or 'right'
    species: str
    chrom: str
    start: int
    end: int
    genes: list[GeneHit]  # ordered by genomic start


# ============================================================================
# Parse
# ============================================================================
def parse_miniprot_gff(path: Path, species: str) -> list[GeneHit]:
    """Parse a whole-genome miniprot GFF, returning all mRNA hits."""
    hits = []
    if not path.exists():
        return hits
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue

            # Strip PanSN prefix if present (species#hap#chrom -> chrom)
            chrom = f[0].split("#", 2)[-1] if "#" in f[0] else f[0]

            # Parse Target=PROTEIN_ID from attributes
            family = "unknown"
            score = 0.0
            for kv in f[8].split(";"):
                if kv.startswith("Target="):
                    family = kv.split("=", 1)[1].split(" ")[0]
                    break
            try:
                score = float(f[5])
            except ValueError:
                pass

            hits.append(GeneHit(
                species=species,
                chrom=chrom,
                start=int(f[3]) - 1,
                end=int(f[4]),
                strand=f[6],
                family=family,
                score=score,
            ))
    return hits


def parse_manifest(path: Path) -> list[dict]:
    species = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 5:
                continue
            species.append({
                "species_id": f[0],
                "genus": f[1],
                "species": f[2],
                "prefix": f[3],
                "fasta": f[4],
                "tier": f[5] if len(f) > 5 else "",
            })
    return species


def parse_bp_bed(path: Path) -> list[dict]:
    """Parse breakpoint BED produced by STEP_05."""
    bps = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            region = f[0]  # e.g. Cgar#1#LG01
            # Strip PanSN
            parts = region.split("#", 2)
            if len(parts) == 3:
                species = parts[0]
                chrom = parts[2]
            else:
                species = "unknown"
                chrom = region
            bps.append({
                "bp_id": f[3],
                "species": species,
                "chrom": chrom,
                "start": int(f[1]),
                "end": int(f[2]),
            })
    return bps


# ============================================================================
# Build indexed gene hits per (species, chrom)
# ============================================================================
def index_by_chrom(hits: list[GeneHit]) -> dict[tuple[str, str], list[GeneHit]]:
    idx: dict[tuple[str, str], list[GeneHit]] = defaultdict(list)
    for h in hits:
        idx[(h.species, h.chrom)].append(h)
    for k in idx:
        idx[k].sort(key=lambda g: g.start)
    return dict(idx)


def extract_flanks(
    bps: list[dict],
    hits_by_chrom: dict[tuple[str, str], list[GeneHit]],
    flank_kb: int,
) -> list[Flank]:
    """
    For each breakpoint, extract left and right flank gene orderings.
    Breakpoint midpoint from BED is (start+end)/2.
    Left flank: [midpoint - flank_kb, midpoint)
    Right flank: [midpoint, midpoint + flank_kb)
    """
    flanks = []
    flank_bp = flank_kb * 1000
    for bp in bps:
        midpoint = (bp["start"] + bp["end"]) // 2
        species = bp["species"]
        chrom = bp["chrom"]
        chrom_hits = hits_by_chrom.get((species, chrom), [])

        left_hits = [h for h in chrom_hits
                     if h.start >= (midpoint - flank_bp) and h.end <= midpoint]
        right_hits = [h for h in chrom_hits
                      if h.start >= midpoint and h.end <= (midpoint + flank_bp)]

        flanks.append(Flank(
            bp_id=bp["bp_id"], side="left",
            species=species, chrom=chrom,
            start=max(0, midpoint - flank_bp), end=midpoint,
            genes=left_hits,
        ))
        flanks.append(Flank(
            bp_id=bp["bp_id"], side="right",
            species=species, chrom=chrom,
            start=midpoint, end=midpoint + flank_bp,
            genes=right_hits,
        ))
    return flanks


# ============================================================================
# REFINE: hamburger-style co-linearity test
# ============================================================================
def find_colinear_cluster_in_chrom(
    query_families: list[str],
    target_chrom_hits: list[GeneHit],
    min_families: int,
    max_gap_genes: int,
    allow_reversed: bool = True,
) -> dict | None:
    """
    Hamburger-style check: does the ordered list of query_families appear as
    a co-linear cluster on the target chromosome?

    Returns the best matching cluster info, or None if no match.

    The target chromosome's gene hits are already sorted by position.
    We walk along target hits and find the longest subsequence of query_families
    that appears in (near-)linear order with bounded gaps.
    """
    if not query_families or not target_chrom_hits:
        return None

    target_families = [g.family for g in target_chrom_hits]

    best = None
    for forward in [True, False] if allow_reversed else [True]:
        qf = list(query_families) if forward else list(reversed(query_families))

        # For each starting position in target, try to match qf against target
        # walking forward with bounded gaps.
        for start_idx in range(len(target_chrom_hits)):
            matched_indices = []
            q_ptr = 0
            t_ptr = start_idx
            gap_since_last_match = 0

            while t_ptr < len(target_chrom_hits) and q_ptr < len(qf):
                if target_families[t_ptr] == qf[q_ptr]:
                    matched_indices.append(t_ptr)
                    q_ptr += 1
                    gap_since_last_match = 0
                else:
                    if matched_indices:
                        gap_since_last_match += 1
                        if gap_since_last_match > max_gap_genes:
                            break
                t_ptr += 1

            if len(matched_indices) >= min_families:
                span = (target_chrom_hits[matched_indices[0]].start,
                        target_chrom_hits[matched_indices[-1]].end)
                record = {
                    "n_matched": len(matched_indices),
                    "fraction": len(matched_indices) / len(query_families),
                    "target_species": target_chrom_hits[0].species,
                    "target_chrom": target_chrom_hits[0].chrom,
                    "target_start": span[0],
                    "target_end": span[1],
                    "orientation": "fwd" if forward else "rev",
                    "matched_families": [qf[i] for i in range(q_ptr)],
                }
                if best is None or record["n_matched"] > best["n_matched"]:
                    best = record

    return best


def refine_flank(
    flank: Flank,
    all_hits_by_chrom: dict[tuple[str, str], list[GeneHit]],
    min_families: int,
    max_gap_genes: int,
) -> dict:
    """
    For this flank, check co-linearity against every other species'
    chromosomes. Return the best hit per other species.
    """
    query_families = [g.family for g in flank.genes if g.family != "unknown"]

    if len(query_families) < min_families:
        return {
            "coherence_class": "TOO_FEW_GENES",
            "query_n_genes": len(flank.genes),
            "query_families": query_families,
            "matches": {},
            "best_match_species": None,
            "best_match_score": 0,
        }

    matches = {}
    for (other_sp, other_chrom), chrom_hits in all_hits_by_chrom.items():
        if other_sp == flank.species:
            continue
        best = find_colinear_cluster_in_chrom(
            query_families, chrom_hits,
            min_families=min_families,
            max_gap_genes=max_gap_genes,
        )
        if best is not None:
            key = other_sp
            if key not in matches or best["n_matched"] > matches[key]["n_matched"]:
                matches[key] = best

    if not matches:
        coh = "LINEAGE_SPECIFIC"
        best_sp = None
        best_score = 0
    else:
        best_sp, best_rec = max(matches.items(), key=lambda kv: kv[1]["n_matched"])
        best_score = best_rec["n_matched"]
        if best_rec["fraction"] >= 0.6:
            coh = "STRONG_ORTHOLOGY"
        elif best_rec["fraction"] >= 0.3:
            coh = "PARTIAL_ORTHOLOGY"
        else:
            coh = "WEAK_ORTHOLOGY"

    return {
        "coherence_class": coh,
        "query_n_genes": len(flank.genes),
        "query_families": query_families,
        "matches": matches,
        "best_match_species": best_sp,
        "best_match_score": best_score,
    }


# ============================================================================
# Per-breakpoint confidence (combine both flanks)
# ============================================================================
def classify_breakpoint_coherence(left: dict, right: dict) -> dict:
    left_class = left["coherence_class"]
    right_class = right["coherence_class"]
    strong = {"STRONG_ORTHOLOGY", "PARTIAL_ORTHOLOGY"}

    if left_class in strong and right_class in strong:
        conf = "HIGH"
        note = "both flanks share orthologous gene order with other species"
    elif (left_class in strong) ^ (right_class in strong):
        conf = "NOVELTY"
        note = ("asymmetric: one flank has preserved gene order, other is "
                "lineage-specific or gene-poor")
    elif left_class == right_class == "LINEAGE_SPECIFIC":
        conf = "LOW"
        note = "both flanks gene-rich but no ortholog match — lineage-specific?"
    elif left_class == right_class == "TOO_FEW_GENES":
        conf = "UNINFORMATIVE"
        note = "both flanks too gene-poor — likely intergenic/repeat desert"
    else:
        conf = "MEDIUM"
        note = "mixed evidence across flanks"

    return {
        "coherence_confidence": conf,
        "coherence_note": note,
        "left_class": left_class,
        "right_class": right_class,
        "left_best_species": left.get("best_match_species", ""),
        "right_best_species": right.get("best_match_species", ""),
        "left_n_matched": left.get("best_match_score", 0),
        "right_n_matched": right.get("best_match_score", 0),
        "left_n_genes": left["query_n_genes"],
        "right_n_genes": right["query_n_genes"],
    }


# ============================================================================
# Main
# ============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bp-bed", required=True, type=Path)
    ap.add_argument("--wg-miniprot-dir", required=True, type=Path,
                    help="Directory with <species_id>.miniprot.gff from STEP_09c1")
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--flank-kb", type=int, default=100,
                    help="Flank size in kb on each side of breakpoint (default 100)")
    ap.add_argument("--min-families", type=int, default=3,
                    help="Min protein families required for co-linearity test (default 3)")
    ap.add_argument("--max-gap-genes", type=int, default=10,
                    help="Max intervening non-matching genes (hamburger-style, default 10)")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    # -------- Load data --------
    print(f"[STEP_09c] Loading manifest: {args.manifest}", file=sys.stderr)
    species_list = parse_manifest(args.manifest)

    print(f"[STEP_09c] Loading breakpoints: {args.bp_bed}", file=sys.stderr)
    bps = parse_bp_bed(args.bp_bed)
    print(f"[STEP_09c]   {len(bps)} breakpoints loaded", file=sys.stderr)

    print(f"[STEP_09c] Loading miniprot GFFs from {args.wg_miniprot_dir}",
          file=sys.stderr)
    all_hits = []
    for sp in species_list:
        gff = args.wg_miniprot_dir / f"{sp['species_id']}.miniprot.gff"
        if not gff.exists():
            print(f"  WARNING: {gff} not found, skipping {sp['species_id']}",
                  file=sys.stderr)
            continue
        hits = parse_miniprot_gff(gff, sp["species_id"])
        print(f"  {sp['species_id']}: {len(hits)} mRNA hits", file=sys.stderr)
        all_hits.extend(hits)

    if not all_hits:
        print("[STEP_09c] ERROR: no miniprot hits loaded. "
              "Did you run STEP_09c1?", file=sys.stderr)
        sys.exit(1)

    hits_by_chrom = index_by_chrom(all_hits)

    # -------- Extract flanks --------
    print(f"\n[STEP_09c] Extracting ±{args.flank_kb} kb flanks around each breakpoint",
          file=sys.stderr)
    flanks = extract_flanks(bps, hits_by_chrom, args.flank_kb)
    flank_gene_counts = [len(f.genes) for f in flanks]
    if flank_gene_counts:
        print(f"  flank gene counts: min={min(flank_gene_counts)}, "
              f"median={sorted(flank_gene_counts)[len(flank_gene_counts)//2]}, "
              f"max={max(flank_gene_counts)}", file=sys.stderr)

    # -------- Refine each flank --------
    print(f"\n[STEP_09c] Refining flanks (co-linearity check, "
          f"min_families={args.min_families}, max_gap={args.max_gap_genes})",
          file=sys.stderr)
    flank_results = {}
    for i, flank in enumerate(flanks):
        if (i + 1) % 20 == 0 or i == 0:
            print(f"  flank {i+1}/{len(flanks)}", file=sys.stderr)
        result = refine_flank(flank, hits_by_chrom,
                              args.min_families, args.max_gap_genes)
        flank_results[(flank.bp_id, flank.side)] = (flank, result)

    # -------- Write per-flank TSV --------
    flank_tsv = args.out / "flank_coherence.tsv"
    with open(flank_tsv, "w") as fh:
        fh.write("bp_id\tside\tspecies\tchrom\tstart\tend\tn_genes\t"
                 "coherence_class\tbest_match_species\tbest_match_score\t"
                 "query_families\n")
        for (bp_id, side), (flank, res) in sorted(flank_results.items()):
            fh.write(f"{bp_id}\t{side}\t{flank.species}\t{flank.chrom}\t"
                     f"{flank.start}\t{flank.end}\t{res['query_n_genes']}\t"
                     f"{res['coherence_class']}\t"
                     f"{res['best_match_species'] or ''}\t"
                     f"{res['best_match_score']}\t"
                     f"{','.join(res['query_families'][:10])}\n")
    print(f"\n[STEP_09c] Per-flank table: {flank_tsv}", file=sys.stderr)

    # -------- Write per-breakpoint confidence --------
    conf_tsv = args.out / "breakpoint_coherence.tsv"
    conf_counts = defaultdict(int)
    with open(conf_tsv, "w") as fh:
        fh.write("bp_id\tcoherence_confidence\tcoherence_note\t"
                 "left_class\tright_class\t"
                 "left_best_species\tright_best_species\t"
                 "left_n_matched\tright_n_matched\t"
                 "left_n_genes\tright_n_genes\n")
        for bp in bps:
            bp_id = bp["bp_id"]
            left_key = (bp_id, "left")
            right_key = (bp_id, "right")
            if left_key not in flank_results or right_key not in flank_results:
                continue
            _, left_res = flank_results[left_key]
            _, right_res = flank_results[right_key]
            cls = classify_breakpoint_coherence(left_res, right_res)
            conf_counts[cls["coherence_confidence"]] += 1
            fh.write(f"{bp_id}\t{cls['coherence_confidence']}\t{cls['coherence_note']}\t"
                     f"{cls['left_class']}\t{cls['right_class']}\t"
                     f"{cls['left_best_species']}\t{cls['right_best_species']}\t"
                     f"{cls['left_n_matched']}\t{cls['right_n_matched']}\t"
                     f"{cls['left_n_genes']}\t{cls['right_n_genes']}\n")

    print(f"[STEP_09c] Per-breakpoint coherence: {conf_tsv}", file=sys.stderr)
    print(f"\n[STEP_09c] Coherence summary:", file=sys.stderr)
    for conf, n in sorted(conf_counts.items()):
        print(f"  {conf:15s} {n:4d}", file=sys.stderr)

    # -------- Write full match details --------
    match_tsv = args.out / "flank_match_details.tsv"
    with open(match_tsv, "w") as fh:
        fh.write("bp_id\tside\ttarget_species\ttarget_chrom\t"
                 "target_start\ttarget_end\torientation\t"
                 "n_matched\tfraction\tmatched_families\n")
        for (bp_id, side), (flank, res) in sorted(flank_results.items()):
            for sp, m in res["matches"].items():
                fh.write(f"{bp_id}\t{side}\t{m['target_species']}\t{m['target_chrom']}\t"
                         f"{m['target_start']}\t{m['target_end']}\t{m['orientation']}\t"
                         f"{m['n_matched']}\t{m['fraction']:.3f}\t"
                         f"{','.join(m['matched_families'])}\n")
    print(f"[STEP_09c] Match details: {match_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
