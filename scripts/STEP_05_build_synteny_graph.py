#!/usr/bin/env python3
"""
STEP_05_build_synteny_graph.py

STAGE B — Synteny graph construction from all-vs-all wfmash PAF.

Builds an undirected multi-species synteny graph where:
    - NODES = genomic intervals bounded by breakpoints (per species, per chromosome)
    - EDGES = homology blocks from wfmash scaffold PAF, with weight = block length
              and an orientation flag (same-strand or inverted)

Produces:
    - nodes.tsv       : interval table (species, chrom, start, end, node_id)
    - edges.tsv       : homology edges (node_u, node_v, length, identity, orientation)
    - breakpoints.bed : inferred breakpoints (one row per species chrom boundary)
    - graph.gexf      : full graph for visualization in Gephi / Cytoscape
    - graph.pickle    : networkx graph for downstream analysis

The graph structure itself encodes phylogenetic signal: a breakpoint shared by
multiple species appears as a shared structural feature in the graph topology.
No pre-specified tree required.

Usage:
    python3 STEP_05_build_synteny_graph.py \
        --paf results/01_wfmash/catfish_clarias.scaffolds.paf \
        --out results/05_synteny_graph/
"""
from __future__ import annotations

import argparse
import pickle
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import networkx as nx


# ----------------------------------------------------------------------------
# Parse PAF
# ----------------------------------------------------------------------------
@dataclass
class Block:
    qsp: str
    qchrom: str
    qstart: int
    qend: int
    strand: str
    tsp: str
    tchrom: str
    tstart: int
    tend: int
    aln_len: int
    matches: int
    identity: float = field(init=False)

    def __post_init__(self):
        self.identity = self.matches / self.aln_len if self.aln_len else 0.0


def parse_paf(path: Path) -> list[Block]:
    blocks = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            qname_parts = f[0].split("#", 2)
            tname_parts = f[5].split("#", 2)
            if len(qname_parts) < 3 or len(tname_parts) < 3:
                continue
            blocks.append(Block(
                qsp=qname_parts[0], qchrom=qname_parts[2],
                qstart=int(f[2]), qend=int(f[3]), strand=f[4],
                tsp=tname_parts[0], tchrom=tname_parts[2],
                tstart=int(f[7]), tend=int(f[8]),
                aln_len=int(f[10]), matches=int(f[9]),
            ))
    return blocks


# ----------------------------------------------------------------------------
# Identify breakpoints per species chromosome
# ----------------------------------------------------------------------------
def identify_breakpoints(
    blocks: list[Block],
    min_block_bp: int = 200_000,
    merge_within_bp: int = 500_000,
) -> dict[tuple[str, str], list[int]]:
    """
    For each (species, chrom), collect all block boundaries that correspond to
    a change in homology — either end of a block mapping to a different target
    chrom, or a strand flip, or a non-collinear jump.

    Returns: {(species, chrom): [sorted breakpoint positions]}
    """
    blocks = [b for b in blocks if (b.qend - b.qstart) >= min_block_bp]

    # For each (species, chrom), collect endpoints where homology changes
    breakpoints: dict[tuple[str, str], set[int]] = defaultdict(set)

    # Group blocks by query (species, chrom)
    by_qchrom: dict[tuple[str, str], list[Block]] = defaultdict(list)
    for b in blocks:
        by_qchrom[(b.qsp, b.qchrom)].append(b)
        # Symmetrically add the target side
        by_qchrom[(b.tsp, b.tchrom)].append(Block(
            qsp=b.tsp, qchrom=b.tchrom, qstart=b.tstart, qend=b.tend,
            strand=b.strand, tsp=b.qsp, tchrom=b.qchrom,
            tstart=b.qstart, tend=b.qend,
            aln_len=b.aln_len, matches=b.matches,
        ))

    for (sp, chrom), blist in by_qchrom.items():
        blist.sort(key=lambda b: b.qstart)
        for i, b in enumerate(blist):
            if i == 0:
                continue
            prev = blist[i - 1]
            # Breakpoint = boundary where homology partner changes
            if (b.tsp, b.tchrom) != (prev.tsp, prev.tchrom) or b.strand != prev.strand:
                # Place breakpoint in the gap between prev.qend and b.qstart
                bp = (prev.qend + b.qstart) // 2
                breakpoints[(sp, chrom)].add(bp)

    # Merge breakpoints within merge_within_bp
    merged: dict[tuple[str, str], list[int]] = {}
    for key, bps in breakpoints.items():
        sorted_bps = sorted(bps)
        if not sorted_bps:
            merged[key] = []
            continue
        out = [sorted_bps[0]]
        for bp in sorted_bps[1:]:
            if bp - out[-1] > merge_within_bp:
                out.append(bp)
            else:
                # Merge: average them
                out[-1] = (out[-1] + bp) // 2
        merged[key] = out

    return merged


# ----------------------------------------------------------------------------
# Build nodes = intervals bounded by breakpoints
# ----------------------------------------------------------------------------
def build_nodes(
    breakpoints: dict[tuple[str, str], list[int]],
    chrom_lengths: dict[tuple[str, str], int],
) -> dict[str, dict]:
    """
    Nodes = intervals on each chromosome, bounded by breakpoints plus chrom
    start/end. Node id format: {species}_{chrom}_{seg_idx}
    """
    nodes: dict[str, dict] = {}

    for (sp, chrom), chrom_len in chrom_lengths.items():
        bps = breakpoints.get((sp, chrom), [])
        boundaries = [0] + sorted(bps) + [chrom_len]
        # Dedupe
        boundaries = sorted(set(boundaries))

        for seg_idx in range(len(boundaries) - 1):
            start = boundaries[seg_idx]
            end = boundaries[seg_idx + 1]
            node_id = f"{sp}_{chrom}_seg{seg_idx:03d}"
            nodes[node_id] = {
                "species": sp,
                "chrom": chrom,
                "start": start,
                "end": end,
                "length": end - start,
                "seg_idx": seg_idx,
            }

    return nodes


def get_chrom_lengths(blocks: list[Block]) -> dict[tuple[str, str], int]:
    """Infer chromosome lengths from the maximum coordinate in PAF."""
    lengths: dict[tuple[str, str], int] = defaultdict(int)
    for b in blocks:
        lengths[(b.qsp, b.qchrom)] = max(lengths[(b.qsp, b.qchrom)], b.qend)
        lengths[(b.tsp, b.tchrom)] = max(lengths[(b.tsp, b.tchrom)], b.tend)
    return dict(lengths)


# ----------------------------------------------------------------------------
# Build edges = homology relationships between nodes
# ----------------------------------------------------------------------------
def find_node_for_position(nodes: dict[str, dict], sp: str, chrom: str, pos: int) -> str | None:
    """Binary-search style — just linear since n_nodes per chrom is small."""
    for node_id, n in nodes.items():
        if n["species"] == sp and n["chrom"] == chrom and n["start"] <= pos < n["end"]:
            return node_id
    return None


def build_edges(
    blocks: list[Block],
    nodes: dict[str, dict],
) -> list[dict]:
    """
    Each block becomes one or more edges between nodes. If a block spans a
    breakpoint (shouldn't happen with proper breakpoint calling but belt-and-
    suspenders), it would generate multiple edges.
    """
    edges = []
    for b in blocks:
        # Find the node containing the midpoint of query and target
        q_mid = (b.qstart + b.qend) // 2
        t_mid = (b.tstart + b.tend) // 2
        q_node = find_node_for_position(nodes, b.qsp, b.qchrom, q_mid)
        t_node = find_node_for_position(nodes, b.tsp, b.tchrom, t_mid)

        if q_node is None or t_node is None:
            continue

        edges.append({
            "u": q_node,
            "v": t_node,
            "length": b.qend - b.qstart,
            "identity": b.identity,
            "strand": b.strand,
            "qstart": b.qstart,
            "qend": b.qend,
            "tstart": b.tstart,
            "tend": b.tend,
        })

    return edges


# ----------------------------------------------------------------------------
# Build networkx graph
# ----------------------------------------------------------------------------
def build_graph(nodes: dict[str, dict], edges: list[dict]) -> nx.MultiGraph:
    G = nx.MultiGraph()
    for node_id, attrs in nodes.items():
        G.add_node(node_id, **attrs)
    for e in edges:
        G.add_edge(e["u"], e["v"],
                   length=e["length"],
                   identity=e["identity"],
                   strand=e["strand"])
    return G


# ----------------------------------------------------------------------------
# Analysis: find structurally shared breakpoints
# ----------------------------------------------------------------------------
def summarize_graph(G: nx.MultiGraph, nodes: dict[str, dict]) -> dict:
    """
    For each node, compute:
      - degree (number of homologous partners)
      - set of species it connects to
      - whether it's a "shared segment" (connects to most species) or
        "lineage-specific" (few connections)
    """
    summary = {}
    all_species = set(n["species"] for n in nodes.values())

    for node_id, attrs in nodes.items():
        partners = set()
        for _, neighbor, _ in G.edges(node_id, keys=True):
            partners.add(nodes[neighbor]["species"])

        summary[node_id] = {
            **attrs,
            "degree": G.degree(node_id),
            "n_species_connected": len(partners),
            "species_set": ",".join(sorted(partners)),
            "is_universal": len(partners) == len(all_species) - 1,
        }

    return summary


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-block-bp", type=int, default=200_000)
    ap.add_argument("--merge-within-bp", type=int, default=500_000,
                    help="Merge breakpoints within this distance")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    print(f"[STEP_05] Loading PAF: {args.paf}", file=sys.stderr)
    blocks = parse_paf(args.paf)
    print(f"[STEP_05] Loaded {len(blocks)} blocks", file=sys.stderr)

    print(f"[STEP_05] Identifying breakpoints...", file=sys.stderr)
    breakpoints = identify_breakpoints(blocks, args.min_block_bp, args.merge_within_bp)
    n_bp = sum(len(v) for v in breakpoints.values())
    print(f"[STEP_05] Found {n_bp} breakpoints across {len(breakpoints)} chromosomes", file=sys.stderr)

    # Write breakpoints BED
    bp_bed = args.out / "breakpoints.bed"
    with open(bp_bed, "w") as fh:
        for (sp, chrom), bps in sorted(breakpoints.items()):
            for i, bp in enumerate(bps):
                fh.write(f"{sp}#1#{chrom}\t{max(0, bp - 100_000)}\t{bp + 100_000}\tbp_{sp}_{chrom}_{i:03d}\n")
    print(f"[STEP_05] Breakpoint BED: {bp_bed}", file=sys.stderr)

    print(f"[STEP_05] Building graph nodes...", file=sys.stderr)
    chrom_lengths = get_chrom_lengths(blocks)
    nodes = build_nodes(breakpoints, chrom_lengths)
    print(f"[STEP_05] Created {len(nodes)} nodes (segments)", file=sys.stderr)

    print(f"[STEP_05] Building graph edges...", file=sys.stderr)
    edges = build_edges(blocks, nodes)
    print(f"[STEP_05] Created {len(edges)} edges", file=sys.stderr)

    print(f"[STEP_05] Building networkx graph...", file=sys.stderr)
    G = build_graph(nodes, edges)

    # Write nodes TSV
    nodes_tsv = args.out / "nodes.tsv"
    with open(nodes_tsv, "w") as fh:
        fh.write("node_id\tspecies\tchrom\tstart\tend\tlength\tseg_idx\n")
        for nid, n in nodes.items():
            fh.write(f"{nid}\t{n['species']}\t{n['chrom']}\t{n['start']}\t{n['end']}\t{n['length']}\t{n['seg_idx']}\n")

    # Write edges TSV
    edges_tsv = args.out / "edges.tsv"
    with open(edges_tsv, "w") as fh:
        fh.write("u\tv\tlength\tidentity\tstrand\tqstart\tqend\ttstart\ttend\n")
        for e in edges:
            fh.write(f"{e['u']}\t{e['v']}\t{e['length']}\t{e['identity']:.4f}\t{e['strand']}\t"
                     f"{e['qstart']}\t{e['qend']}\t{e['tstart']}\t{e['tend']}\n")

    # Graph summary
    summary = summarize_graph(G, nodes)
    summary_tsv = args.out / "node_summary.tsv"
    with open(summary_tsv, "w") as fh:
        fh.write("node_id\tspecies\tchrom\tstart\tend\tlength\tdegree\tn_species_connected\tspecies_set\tis_universal\n")
        for nid, s in summary.items():
            fh.write(f"{nid}\t{s['species']}\t{s['chrom']}\t{s['start']}\t{s['end']}\t"
                     f"{s['length']}\t{s['degree']}\t{s['n_species_connected']}\t"
                     f"{s['species_set']}\t{s['is_universal']}\n")

    # Write graph files
    graph_gexf = args.out / "graph.gexf"
    # networkx GEXF doesn't love multigraphs with complex attrs; convert to simple
    G_simple = nx.Graph()
    G_simple.add_nodes_from((n, {k: v for k, v in d.items()}) for n, d in G.nodes(data=True))
    for u, v, data in G.edges(data=True):
        if G_simple.has_edge(u, v):
            G_simple[u][v]["total_length"] += data["length"]
            G_simple[u][v]["n_blocks"] += 1
        else:
            G_simple.add_edge(u, v,
                              total_length=data["length"],
                              n_blocks=1,
                              mean_identity=data["identity"],
                              strand=data["strand"])
    nx.write_gexf(G_simple, graph_gexf)
    print(f"[STEP_05] GEXF (for Gephi/Cytoscape): {graph_gexf}", file=sys.stderr)

    # Pickle for downstream
    with open(args.out / "graph.pickle", "wb") as fh:
        pickle.dump({"graph": G, "nodes": nodes, "edges": edges,
                     "breakpoints": breakpoints}, fh)

    print(f"\n[STEP_05] Graph summary:", file=sys.stderr)
    print(f"  Nodes:                  {G.number_of_nodes()}", file=sys.stderr)
    print(f"  Edges (homology links): {G.number_of_edges()}", file=sys.stderr)
    print(f"  Connected components:   {nx.number_connected_components(G)}", file=sys.stderr)
    universal = sum(1 for s in summary.values() if s["is_universal"])
    print(f"  'Universal' segments (connect to all species): {universal}", file=sys.stderr)
    singletons = sum(1 for s in summary.values() if s["n_species_connected"] == 0)
    print(f"  Orphan segments (no homology): {singletons}", file=sys.stderr)


if __name__ == "__main__":
    main()
