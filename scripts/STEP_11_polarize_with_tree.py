#!/usr/bin/env python3
"""
STEP_11_polarize_with_tree.py

Use a species tree (Newick format) to polarize breakpoints detected in Stage B.

Given:
    - Your BUSCO supermatrix species tree (Newick)
    - Per-breakpoint flank coherence (STEP_09c flank_match_details.tsv)

For each breakpoint:
    1. Look at which species share the flank gene order at each side
    2. Map those groupings onto the species tree
    3. Use parsimony to infer on which branch the rearrangement most likely
       occurred (the "polarized direction")

This gives you the actual evolutionary direction of each rearrangement
("fusion in Cgar lineage" vs "fission in Cmac lineage") IF the tree topology
supports a unique polarization. When the pattern is ambiguous, the script
reports that honestly rather than forcing a direction.

Requires: ete3 (pip install ete3)
Usage:
    python3 STEP_11_polarize_with_tree.py \\
        --tree config/species_tree.nwk \\
        --flank-details results/09c_flank_coherence/flank_match_details.tsv \\
        --bp-coherence results/09c_flank_coherence/breakpoint_coherence.tsv \\
        --out results/11_polarized_breakpoints/
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def load_tree_newick(path: Path) -> "ete3.Tree":
    try:
        from ete3 import Tree
    except ImportError:
        print("ERROR: ete3 not installed. Run: pip install ete3", file=sys.stderr)
        sys.exit(1)
    return Tree(str(path), format=1)


def get_mrca_species(tree, species_set: set[str]) -> str | None:
    """Return the MRCA node name for a set of species leaves."""
    try:
        mrca = tree.get_common_ancestor(list(species_set))
        if mrca.name:
            return mrca.name
        leaves = sorted([l.name for l in mrca.get_leaves()])
        return "MRCA_of_" + "+".join(leaves)
    except Exception:
        return None


def polarize_breakpoint(
    bp_id: str,
    left_matches: dict[str, set[str]],   # left_flank species: set of target species sharing flank
    right_matches: dict[str, set[str]],
    focal_species: str,
    tree,
) -> dict:
    """
    Parsimony-style polarization.

    For each flank, collect the species that share the flank's gene order
    (the 'ingroup' for that flank). The branch on the species tree where the
    rearrangement occurred is the one separating the sharers from non-sharers.
    """
    all_leaves = set(l.name for l in tree.get_leaves())

    left_sharers = set(left_matches.get(focal_species, set())) | {focal_species}
    right_sharers = set(right_matches.get(focal_species, set())) | {focal_species}

    left_sharers &= all_leaves
    right_sharers &= all_leaves

    L_mrca = get_mrca_species(tree, left_sharers) if len(left_sharers) >= 2 else None
    R_mrca = get_mrca_species(tree, right_sharers) if len(right_sharers) >= 2 else None

    # Polarization logic:
    if left_sharers == right_sharers:
        polarization = "SHARED_STATE"
        note = "both flanks share the same species group — rearrangement fixed in this group"
        branch = L_mrca
    elif left_sharers.issubset(right_sharers) and left_sharers != right_sharers:
        polarization = "LEFT_DERIVED"
        note = f"right flank is ancestral (broader sharing); left flank rearranged on branch to {focal_species}"
        branch = focal_species
    elif right_sharers.issubset(left_sharers) and left_sharers != right_sharers:
        polarization = "RIGHT_DERIVED"
        note = f"left flank is ancestral (broader sharing); right flank rearranged on branch to {focal_species}"
        branch = focal_species
    elif not left_sharers.intersection(right_sharers) - {focal_species}:
        polarization = "AMBIGUOUS"
        note = "flanks share with disjoint species sets — likely independent events or complex rearrangement"
        branch = None
    else:
        polarization = "PARTIAL_OVERLAP"
        note = f"flanks share with overlapping but not nested species sets"
        branch = None

    return {
        "bp_id": bp_id,
        "focal_species": focal_species,
        "left_sharers": ",".join(sorted(left_sharers)),
        "right_sharers": ",".join(sorted(right_sharers)),
        "left_mrca": L_mrca or "",
        "right_mrca": R_mrca or "",
        "polarization": polarization,
        "inferred_branch": branch or "",
        "note": note,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True, type=Path, help="Newick species tree")
    ap.add_argument("--flank-details", required=True, type=Path,
                    help="STEP_09c flank_match_details.tsv")
    ap.add_argument("--bp-coherence", required=True, type=Path,
                    help="STEP_09c breakpoint_coherence.tsv (for focal species lookup)")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    print(f"[STEP_11] Loading species tree: {args.tree}", file=sys.stderr)
    tree = load_tree_newick(args.tree)
    leaves = [l.name for l in tree.get_leaves()]
    print(f"[STEP_11]   {len(leaves)} leaves: {', '.join(leaves)}", file=sys.stderr)

    # Load flank match details
    # Columns: bp_id, side, target_species, target_chrom, ...
    flank_matches_by_bp: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"left": set(), "right": set()}
    )
    with open(args.flank_details) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            bp = row["bp_id"]
            side = row["side"]
            tsp = row["target_species"]
            flank_matches_by_bp[bp][side].add(tsp)

    # Focal species from bp_coherence file
    focal_by_bp: dict[str, str] = {}
    with open(args.bp_coherence) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # bp_id format "bp_<species>_<chrom>_<idx>" — the species part is position [1]
            bp = row["bp_id"]
            parts = bp.split("_")
            focal_by_bp[bp] = parts[1] if len(parts) > 1 else ""

    # Polarize
    print(f"[STEP_11] Polarizing {len(focal_by_bp)} breakpoints", file=sys.stderr)
    results = []
    for bp_id, focal in focal_by_bp.items():
        if not focal:
            continue
        lm = {focal: flank_matches_by_bp[bp_id]["left"]}
        rm = {focal: flank_matches_by_bp[bp_id]["right"]}
        r = polarize_breakpoint(bp_id, lm, rm, focal, tree)
        results.append(r)

    # Write output
    out_tsv = args.out / "polarized_breakpoints.tsv"
    with open(out_tsv, "w") as fh:
        fields = ["bp_id", "focal_species", "left_sharers", "right_sharers",
                  "left_mrca", "right_mrca", "polarization", "inferred_branch", "note"]
        fh.write("\t".join(fields) + "\n")
        for r in results:
            fh.write("\t".join(str(r.get(f, "")) for f in fields) + "\n")

    # Summary
    counts = defaultdict(int)
    for r in results:
        counts[r["polarization"]] += 1
    print(f"\n[STEP_11] Polarization summary:", file=sys.stderr)
    for k, v in sorted(counts.items()):
        print(f"  {k:20s} {v:4d}", file=sys.stderr)
    print(f"\n[STEP_11] Output: {out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
