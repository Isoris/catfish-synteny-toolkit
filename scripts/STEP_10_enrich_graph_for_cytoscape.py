#!/usr/bin/env python3
"""
STEP_10_enrich_graph_for_cytoscape.py

Consume all downstream analysis outputs and attach them as node/edge
attributes on the synteny graph from STEP_05. Produces an enriched GEXF
file ready for Cytoscape filtering and coloring.

Attributes attached per NODE (segment between breakpoints):
  Positive evidence:
    support_kmer_median_identity  : median ANI of edges touching this node
    support_protein_coherence     : STRONG / PARTIAL / WEAK / LINEAGE_SPECIFIC
    support_orthology_score       : 0-3 (count of agreeing evidence layers)
    n_species_connected           : preserved from STEP_05

  Negative evidence:
    te_fraction                   : fraction of node covered by TEs (if available)
    satellite_bp                  : TRF satellite bp inside this segment
    palindrome_count              : IRF inverted repeat count
    sd_count                      : minimap2 self-alignment hit count
    is_assembly_suspect           : True if te+satellite > 50% of segment

Attributes per EDGE (homology link):
  identity, length, strand        : preserved from STEP_05
  evidence_agreement              : count of independent evidence types
                                    that support this edge (0-3)

Usage:
    python3 STEP_10_enrich_graph_for_cytoscape.py \\
        --graph-pickle results/05_synteny_graph/graph.pickle \\
        --bp-annotation results/07_bp_annotation/breakpoint_annotation_summary.tsv \\
        --flank-coherence results/09c_flank_coherence/flank_coherence.tsv \\
        --dual-evidence results/09_dual_evidence/breakpoint_confidence.tsv \\
        --out results/10_enriched_graph/
"""
from __future__ import annotations

import argparse
import csv
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import networkx as nx


# ----------------------------------------------------------------------------
# Load attribute sources
# ----------------------------------------------------------------------------
def load_tsv(path: Path, key_field: str) -> dict[str, dict]:
    """Parse TSV with header; key on the given field."""
    if path is None or not path.exists():
        return {}
    out = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if key_field in row:
                out[row[key_field]] = row
    return out


# ----------------------------------------------------------------------------
# Attach node attributes by finding breakpoints within each node
# ----------------------------------------------------------------------------
def enrich_nodes(
    G: nx.MultiGraph,
    nodes_meta: dict,
    bp_annot: dict[str, dict],
    flank_coh: dict[tuple[str, str], dict],
    bp_conf: dict[str, dict],
) -> None:
    """
    For each graph node (a segment), find the breakpoints that bound it,
    pull in the annotation for those breakpoints, compute aggregate scores.
    """
    for node_id, attrs in G.nodes(data=True):
        species = attrs.get("species", "")
        chrom = attrs.get("chrom", "")
        start = attrs.get("start", 0)
        end = attrs.get("end", 0)

        # Find breakpoints that touch this segment (at either boundary)
        bp_ids_at_boundaries = []
        for bp_id, row in bp_annot.items():
            bp_id_species = bp_id.split("_")[1] if "_" in bp_id else ""
            if bp_id_species == species:
                # bp_id format: bp_<species>_<chrom>_<idx>
                # region info might be in a separate field; use 'region' col
                region = row.get("region", "")
                if region:
                    # "Cgar#1#LG01:950000-1050000"
                    try:
                        _, coords = region.split(":", 1)
                        bp_start, bp_end = coords.split("-")
                        bp_mid = (int(bp_start) + int(bp_end)) // 2
                        if abs(bp_mid - start) < 100_000 or abs(bp_mid - end) < 100_000:
                            bp_ids_at_boundaries.append(bp_id)
                    except Exception:
                        pass

        # Aggregate annotation from boundaries
        te_pct_vals = []
        sat_bp_total = 0
        palindrome_total = 0
        sd_total = 0
        for bp_id in bp_ids_at_boundaries:
            row = bp_annot.get(bp_id, {})
            try:
                te_pct_vals.append(float(row.get("te_pct", 0) or 0))
            except ValueError:
                pass
            try:
                sat_bp_total += int(row.get("trf_satellite_bp", 0) or 0)
                palindrome_total += int(row.get("irf_n_palindromes", 0) or 0)
                sd_total += int(row.get("selfaln_n_hits", 0) or 0)
            except ValueError:
                pass
        te_fraction = (sum(te_pct_vals) / len(te_pct_vals) / 100.0) if te_pct_vals else 0.0

        # Aggregate flank coherence from STEP_09c (by bp_id + side)
        orthology_classes = []
        for bp_id in bp_ids_at_boundaries:
            for side in ("left", "right"):
                row = flank_coh.get((bp_id, side)) if isinstance(flank_coh, dict) else None
                if row:
                    orthology_classes.append(row.get("coherence_class", ""))

        protein_coherence = "UNKNOWN"
        if "STRONG_ORTHOLOGY" in orthology_classes:
            protein_coherence = "STRONG"
        elif "PARTIAL_ORTHOLOGY" in orthology_classes:
            protein_coherence = "PARTIAL"
        elif "WEAK_ORTHOLOGY" in orthology_classes:
            protein_coherence = "WEAK"
        elif "LINEAGE_SPECIFIC" in orthology_classes:
            protein_coherence = "LINEAGE_SPECIFIC"
        elif "TOO_FEW_GENES" in orthology_classes:
            protein_coherence = "GENE_POOR"

        # K-mer identity from edges
        ident_vals = []
        for _, _, d in G.edges(node_id, data=True):
            if "identity" in d:
                try:
                    ident_vals.append(float(d["identity"]))
                except (TypeError, ValueError):
                    pass
        kmer_median = sum(ident_vals) / len(ident_vals) if ident_vals else 0.0

        # Dual-evidence confidence from STEP_09b
        dual_conf = "UNKNOWN"
        for bp_id in bp_ids_at_boundaries:
            row = bp_conf.get(bp_id, {})
            if row:
                dual_conf = row.get("confidence", "UNKNOWN")
                break

        # Compute aggregate orthology score (0-3 layers agreeing)
        score = 0
        if kmer_median >= 0.85:
            score += 1
        if dual_conf == "HIGH":
            score += 1
        if protein_coherence in ("STRONG", "PARTIAL"):
            score += 1

        # Assembly suspect flag
        segment_len = max(1, end - start)
        suspicious_fraction = te_fraction + (sat_bp_total / segment_len)
        is_suspect = suspicious_fraction > 0.5

        # Attach to node
        G.nodes[node_id]["support_kmer_identity"] = round(kmer_median, 4)
        G.nodes[node_id]["support_protein_coherence"] = protein_coherence
        G.nodes[node_id]["support_dual_confidence"] = dual_conf
        G.nodes[node_id]["support_orthology_score"] = score
        G.nodes[node_id]["te_fraction"] = round(te_fraction, 3)
        G.nodes[node_id]["satellite_bp"] = sat_bp_total
        G.nodes[node_id]["palindrome_count"] = palindrome_total
        G.nodes[node_id]["sd_count"] = sd_total
        G.nodes[node_id]["is_assembly_suspect"] = is_suspect
        G.nodes[node_id]["n_bp_at_boundaries"] = len(bp_ids_at_boundaries)


# ----------------------------------------------------------------------------
# Attach edge attributes
# ----------------------------------------------------------------------------
def enrich_edges(G: nx.MultiGraph, nodes_meta: dict) -> None:
    """
    For each edge, compute an 'evidence_agreement' score based on the
    flanking nodes' attributes.
    """
    for u, v, k, d in G.edges(keys=True, data=True):
        u_attrs = G.nodes[u]
        v_attrs = G.nodes[v]

        agreement = 0
        if d.get("identity", 0) >= 0.85:
            agreement += 1
        if (u_attrs.get("support_protein_coherence") in ("STRONG", "PARTIAL")
                and v_attrs.get("support_protein_coherence") in ("STRONG", "PARTIAL")):
            agreement += 1
        if (not u_attrs.get("is_assembly_suspect", False)
                and not v_attrs.get("is_assembly_suspect", False)):
            agreement += 1

        G.edges[u, v, k]["evidence_agreement"] = agreement


# ----------------------------------------------------------------------------
# Write simple-graph GEXF (Cytoscape doesn't love multigraphs)
# ----------------------------------------------------------------------------
def write_simple_gexf(G: nx.MultiGraph, out_path: Path) -> None:
    """Collapse multigraph into simple graph, aggregating parallel edges."""
    S = nx.Graph()
    for n, attrs in G.nodes(data=True):
        # GEXF requires scalar attrs; convert booleans and drop None
        clean = {}
        for k, v in attrs.items():
            if v is None:
                continue
            if isinstance(v, bool):
                clean[k] = "true" if v else "false"
            else:
                clean[k] = v
        S.add_node(n, **clean)

    for u, v, d in G.edges(data=True):
        if S.has_edge(u, v):
            S[u][v]["n_links"] += 1
            S[u][v]["total_length"] += d.get("length", 0)
            S[u][v]["max_agreement"] = max(
                S[u][v]["max_agreement"], d.get("evidence_agreement", 0)
            )
        else:
            S.add_edge(u, v,
                       n_links=1,
                       total_length=int(d.get("length", 0)),
                       mean_identity=float(d.get("identity", 0.0)),
                       strand=str(d.get("strand", "+")),
                       max_agreement=int(d.get("evidence_agreement", 0)))

    nx.write_gexf(S, out_path)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph-pickle", required=True, type=Path)
    ap.add_argument("--bp-annotation", type=Path, default=None,
                    help="STEP_07 breakpoint_annotation_summary.tsv")
    ap.add_argument("--flank-coherence", type=Path, default=None,
                    help="STEP_09c flank_coherence.tsv")
    ap.add_argument("--dual-evidence", type=Path, default=None,
                    help="STEP_09b breakpoint_confidence.tsv")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    # Load graph
    print(f"[STEP_10] Loading graph: {args.graph_pickle}", file=sys.stderr)
    with open(args.graph_pickle, "rb") as fh:
        data = pickle.load(fh)
    G = data["graph"]
    nodes_meta = data["nodes"]
    print(f"[STEP_10]   {G.number_of_nodes()} nodes, {G.number_of_edges()} edges",
          file=sys.stderr)

    # Load attribute sources
    print(f"[STEP_10] Loading annotation sources...", file=sys.stderr)
    bp_annot = load_tsv(args.bp_annotation, key_field="bp_id") if args.bp_annotation else {}
    bp_conf = load_tsv(args.dual_evidence, key_field="bp_id") if args.dual_evidence else {}
    print(f"  bp_annot rows:  {len(bp_annot)}", file=sys.stderr)
    print(f"  bp_conf rows:   {len(bp_conf)}", file=sys.stderr)

    # Flank coherence: two rows per bp (left, right) — key on (bp_id, side)
    flank_coh = {}
    if args.flank_coherence and args.flank_coherence.exists():
        with open(args.flank_coherence) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                flank_coh[(row["bp_id"], row["side"])] = row
    print(f"  flank_coh rows: {len(flank_coh)}", file=sys.stderr)

    # Enrich
    print(f"[STEP_10] Enriching nodes...", file=sys.stderr)
    enrich_nodes(G, nodes_meta, bp_annot, flank_coh, bp_conf)
    print(f"[STEP_10] Enriching edges...", file=sys.stderr)
    enrich_edges(G, nodes_meta)

    # Write outputs
    enriched_gexf = args.out / "graph_enriched.gexf"
    write_simple_gexf(G, enriched_gexf)
    print(f"[STEP_10] Enriched GEXF (for Cytoscape): {enriched_gexf}", file=sys.stderr)

    # Also dump a flat node TSV with all scores for quick grepping
    node_tsv = args.out / "nodes_enriched.tsv"
    with open(node_tsv, "w") as fh:
        fields = ["node_id", "species", "chrom", "start", "end", "length",
                 "n_species_connected",
                 "support_kmer_identity", "support_protein_coherence",
                 "support_dual_confidence", "support_orthology_score",
                 "te_fraction", "satellite_bp", "palindrome_count", "sd_count",
                 "is_assembly_suspect", "n_bp_at_boundaries"]
        fh.write("\t".join(fields) + "\n")
        for node_id, attrs in G.nodes(data=True):
            vals = [str(attrs.get(f, "") if f != "node_id" else node_id) for f in fields]
            fh.write("\t".join(vals) + "\n")
    print(f"[STEP_10] Enriched node TSV: {node_tsv}", file=sys.stderr)

    # Summary
    print(f"\n[STEP_10] Score distributions:", file=sys.stderr)
    score_counts = defaultdict(int)
    for _, attrs in G.nodes(data=True):
        score_counts[attrs.get("support_orthology_score", 0)] += 1
    for s in sorted(score_counts):
        print(f"  orthology_score={s}: {score_counts[s]} nodes", file=sys.stderr)

    suspect_n = sum(1 for _, a in G.nodes(data=True) if a.get("is_assembly_suspect"))
    print(f"  Assembly-suspect nodes: {suspect_n}/{G.number_of_nodes()}", file=sys.stderr)


if __name__ == "__main__":
    main()
