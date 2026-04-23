#!/usr/bin/env python3
"""
STEP_06_plot_synteny_graph.py

STAGE B continued — visualize the synteny graph.

Produces three views:
    1. Linear chromosome ladder with edges (Kuang Figure 1B-style ribbons,
       but built from the graph instead of raw PAF)
    2. Force-directed graph layout (Gephi-style, for topology inspection)
    3. Species-pair adjacency matrix heatmap (shared-segment counts)

Usage:
    python3 STEP_06_plot_synteny_graph.py --graph-dir results/05_synteny_graph/
"""
from __future__ import annotations

import argparse
import pickle
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np


def load_graph(graph_dir: Path) -> dict:
    with open(graph_dir / "graph.pickle", "rb") as fh:
        return pickle.load(fh)


# ----------------------------------------------------------------------------
# Plot 1: Linear synteny ribbons from the graph
# ----------------------------------------------------------------------------
def plot_linear_ribbons(data: dict, out_path: Path, chrom_gap: int = 5_000_000):
    G = data["graph"]
    nodes = data["nodes"]

    # Determine species order and chromosome order per species
    species_list = sorted(set(n["species"] for n in nodes.values()))

    # Chromosome layout: stack species as rows, chroms side-by-side per row
    chrom_positions: dict[tuple[str, str], tuple[int, int]] = {}  # (sp, chrom) -> (x_start, x_end)
    chrom_y: dict[str, float] = {}

    def chrom_sort_key(chrom):
        # extract trailing number if present
        import re
        m = re.search(r"\d+$", chrom)
        return (int(m.group()) if m else 999999, chrom)

    for sp_idx, sp in enumerate(species_list):
        chrom_y[sp] = len(species_list) - sp_idx  # top = first species
        sp_chroms = sorted(
            set(n["chrom"] for n in nodes.values() if n["species"] == sp),
            key=chrom_sort_key,
        )
        x_cursor = 0
        for chrom in sp_chroms:
            chrom_len = max(
                n["end"] for n in nodes.values()
                if n["species"] == sp and n["chrom"] == chrom
            )
            chrom_positions[(sp, chrom)] = (x_cursor, x_cursor + chrom_len)
            x_cursor += chrom_len + chrom_gap

    fig, ax = plt.subplots(figsize=(20, 2 + 1.2 * len(species_list)))

    # Draw chromosome bars
    for (sp, chrom), (xs, xe) in chrom_positions.items():
        y = chrom_y[sp]
        ax.add_patch(mpatches.Rectangle(
            (xs, y - 0.1), xe - xs, 0.2,
            facecolor="lightgrey", edgecolor="black", linewidth=0.3,
        ))
        # Chrom label at midpoint
        ax.text((xs + xe) / 2, y + 0.15, chrom, ha="center", va="bottom",
                fontsize=5, rotation=0)

    # Draw ribbons from edges
    drawn = set()
    for u, v, data_e in G.edges(data=True):
        if (u, v) in drawn or (v, u) in drawn:
            continue
        drawn.add((u, v))

        nu = nodes[u]
        nv = nodes[v]
        if nu["species"] == nv["species"]:
            continue  # skip intra-species

        xu_start, _ = chrom_positions[(nu["species"], nu["chrom"])]
        xv_start, _ = chrom_positions[(nv["species"], nv["chrom"])]

        x_u = xu_start + (nu["start"] + nu["end"]) / 2
        x_v = xv_start + (nv["start"] + nv["end"]) / 2
        y_u = chrom_y[nu["species"]]
        y_v = chrom_y[nv["species"]]

        color = "red" if data_e.get("strand") == "-" else "steelblue"
        ax.plot([x_u, x_v], [y_u, y_v], color=color, alpha=0.1, linewidth=0.4)

    # Species labels on left
    ax.set_yticks(list(chrom_y.values()))
    ax.set_yticklabels(list(chrom_y.keys()), fontsize=9)
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_title("Synteny ribbons across species (from graph)")
    ax.spines[['top', 'right', 'bottom']].set_visible(False)

    # Legend
    ax.plot([], [], color="steelblue", label="Forward")
    ax.plot([], [], color="red", label="Inverted")
    ax.legend(loc="upper right", fontsize=8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


# ----------------------------------------------------------------------------
# Plot 2: Species-pair adjacency matrix (shared segments)
# ----------------------------------------------------------------------------
def plot_pair_heatmap(data: dict, out_path: Path):
    G = data["graph"]
    nodes = data["nodes"]

    species_list = sorted(set(n["species"] for n in nodes.values()))
    n = len(species_list)
    sp_idx = {sp: i for i, sp in enumerate(species_list)}

    mat = np.zeros((n, n), dtype=float)
    for u, v, d in G.edges(data=True):
        s1 = nodes[u]["species"]
        s2 = nodes[v]["species"]
        if s1 != s2:
            mat[sp_idx[s1], sp_idx[s2]] += d["length"]
            if s1 != s2:
                mat[sp_idx[s2], sp_idx[s1]] += d["length"]

    # Convert to Mb
    mat = mat / 1e6

    fig, ax = plt.subplots(figsize=(1 + 0.5 * n, 1 + 0.5 * n))
    im = ax.imshow(mat, cmap="viridis")
    ax.set_xticks(range(n)); ax.set_xticklabels(species_list, rotation=45, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(species_list)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{mat[i,j]:.0f}", ha="center", va="center",
                    color="white" if mat[i, j] < mat.max() / 2 else "black", fontsize=7)
    ax.set_title("Pairwise syntenic sequence (Mb)")
    plt.colorbar(im, ax=ax, label="Shared synteny (Mb)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


# ----------------------------------------------------------------------------
# Plot 3: Breakpoint density per species
# ----------------------------------------------------------------------------
def plot_breakpoint_density(data: dict, out_path: Path):
    breakpoints = data["breakpoints"]

    # Count breakpoints per species
    counts = defaultdict(int)
    for (sp, _), bps in breakpoints.items():
        counts[sp] += len(bps)

    species_list = sorted(counts.keys())
    values = [counts[sp] for sp in species_list]

    fig, ax = plt.subplots(figsize=(max(6, 0.8 * len(species_list)), 4))
    ax.bar(species_list, values, color="steelblue", edgecolor="black")
    ax.set_ylabel("Number of breakpoints")
    ax.set_title("Breakpoints per species (from synteny graph)")
    ax.tick_params(axis='x', rotation=45)
    for i, v in enumerate(values):
        ax.text(i, v + 0.5, str(v), ha="center", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--graph-dir", required=True, type=Path)
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Output dir (default: <graph-dir>/figures)")
    args = ap.parse_args()

    out_dir = args.out_dir or (args.graph_dir / "figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[STEP_06] Loading graph from {args.graph_dir}")
    data = load_graph(args.graph_dir)

    print(f"[STEP_06] Plotting linear ribbons...")
    plot_linear_ribbons(data, out_dir / "synteny_ribbons.pdf")

    print(f"[STEP_06] Plotting pair heatmap...")
    plot_pair_heatmap(data, out_dir / "pair_heatmap.pdf")

    print(f"[STEP_06] Plotting breakpoint density...")
    plot_breakpoint_density(data, out_dir / "breakpoint_density.pdf")

    print(f"[STEP_06] Done. Figures in: {out_dir}")


if __name__ == "__main__":
    main()
