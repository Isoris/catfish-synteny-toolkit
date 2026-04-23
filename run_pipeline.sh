#!/usr/bin/env bash
# ============================================================================
# run_pipeline.sh — master orchestration for catfish-synteny-toolkit
#
# Runs the full pipeline in the correct order. Submits SLURM jobs with
# dependencies. Safe to re-run: skips stages whose outputs already exist
# unless --force is specified.
#
# Usage:
#   bash run_pipeline.sh <tier> [--stage-start N] [--stage-end N] [--force]
#
# Stages (each is one numbered script set):
#   0  = PanSN concat
#   1  = wfmash all-vs-all
#   5  = synteny graph
#   6  = graph figures
#   7  = breakpoint annotation (TRF/IRF/SDs/TE)
#   8  = cross-species satellite stats
#   9  = dual-evidence validation (kmer + per-protein)
#   9c = flank coherence (protein-family order)
#   10 = enriched graph for Cytoscape
#   11 = tree-based polarization
#
# Example — run Stage A only (fast scoping):
#   bash run_pipeline.sh core --stage-end 1
#
# Example — run full pipeline for Clarias tier:
#   bash run_pipeline.sh clarias
#
# Example — redo just annotations (assumes upstream done):
#   bash run_pipeline.sh clarias --stage-start 7 --stage-end 8
# ============================================================================
set -euo pipefail

TIER="${1:?Usage: run_pipeline.sh <tier: core|clarias|all|deep> [options]}"
shift

STAGE_START=0
STAGE_END=11
FORCE=0
TE_LIB=""
PROTEOME=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --stage-start) STAGE_START="$2"; shift 2 ;;
        --stage-end)   STAGE_END="$2"; shift 2 ;;
        --force)       FORCE=1; shift ;;
        --te-lib)      TE_LIB="$2"; shift 2 ;;
        --proteome)    PROTEOME="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

mkdir -p logs
SCRIPTS="scripts"

stage_done() {
    local marker="$1"
    [[ -f "$marker" && $FORCE -eq 0 ]]
}

run_stage() {
    local n="$1"
    local name="$2"
    local cmd="$3"
    local marker="$4"

    if [[ $n -lt $STAGE_START || $n -gt $STAGE_END ]]; then
        echo "[skip]  Stage $n ($name) — outside requested range"
        return 0
    fi
    if stage_done "$marker"; then
        echo "[cache] Stage $n ($name) — marker exists, skip (use --force to rerun)"
        return 0
    fi
    echo ""
    echo "=== Stage $n: $name ==="
    echo "    $cmd"
    eval "$cmd"
    touch "$marker"
}

# ---- STAGE 0b: per-species ScaffoldKit prep (SLURM array) ----
if [[ $STAGE_START -le 0 && $STAGE_END -ge 0 ]] && ! stage_done "results/00_prepared_genomes/.done"; then
    echo ""
    echo "=== Stage 0b: per-species ScaffoldKit preparation (SLURM array) ==="
    N_SP=$(tail -n +2 config/species_manifest.tsv | awk -F'\t' -v pat="^(core|clarias_context|far_outgroup|deep_outgroup)" '$6 ~ pat && $1 != "" && $1 !~ /^#/' | wc -l)
    if (( N_SP > 0 )); then
        JOBID=$(sbatch --parsable --array=0-$((N_SP-1)) $SCRIPTS/STEP_00b_prepare_genomes.sh $TIER)
        echo "    Submitted array $JOBID; waiting..."
        while squeue -j "$JOBID" &>/dev/null; do sleep 60; done
    fi
    mkdir -p results/00_prepared_genomes
    touch results/00_prepared_genomes/.done
fi

# ---- STAGE 0 ----
run_stage 0 "PanSN concat" \
    "bash $SCRIPTS/STEP_00_prepare_panSN.sh $TIER" \
    "results/00_panSN_inputs/catfish_${TIER}.fa.gz.fai"

# ---- STAGE 1 (SLURM — blocks until done) ----
if [[ $STAGE_START -le 1 && $STAGE_END -ge 1 ]] && ! stage_done "results/01_wfmash/catfish_${TIER}.scaffolds.paf"; then
    echo ""
    echo "=== Stage 1: wfmash all-vs-all (submitting SLURM) ==="
    JOBID=$(sbatch --parsable $SCRIPTS/STEP_01_wfmash_allvsall.sh $TIER approx)
    echo "    Submitted job $JOBID; waiting for completion..."
    while squeue -j "$JOBID" &>/dev/null; do sleep 60; done
    [[ -f "results/01_wfmash/catfish_${TIER}.scaffolds.paf" ]] || { echo "STAGE 1 FAILED"; exit 1; }
fi

# ---- STAGE 5: synteny graph ----
run_stage 5 "Synteny graph" \
    "python3 $SCRIPTS/STEP_05_build_synteny_graph.py \
        --paf results/01_wfmash/catfish_${TIER}.scaffolds.paf \
        --out results/05_synteny_graph/" \
    "results/05_synteny_graph/breakpoints.bed"

# ---- STAGE 6: graph figures ----
run_stage 6 "Graph figures" \
    "python3 $SCRIPTS/STEP_06_plot_synteny_graph.py --graph-dir results/05_synteny_graph/" \
    "results/05_synteny_graph/figures/synteny_ribbons.pdf"

echo ""
echo "=============================================================="
echo "Decision checkpoint — look at these files before continuing:"
echo "  results/05_synteny_graph/figures/synteny_ribbons.pdf"
echo "  results/05_synteny_graph/figures/pair_heatmap.pdf"
echo "  results/05_synteny_graph/figures/breakpoint_density.pdf"
echo ""
echo "If the signal is weak, STOP here. Re-focus on the inversion paper."
echo "If signal is clear, continue with --stage-start 7"
echo "=============================================================="

if [[ $STAGE_END -lt 7 ]]; then
    exit 0
fi

# ---- STAGE 7: breakpoint annotation (SLURM) ----
if [[ $STAGE_START -le 7 && $STAGE_END -ge 7 ]] && ! stage_done "results/07_bp_annotation/breakpoint_annotation_summary.tsv"; then
    echo ""
    echo "=== Stage 7: breakpoint annotation (submitting SLURM) ==="
    TE_ARG=""
    [[ -n "$TE_LIB" ]] && TE_ARG="$TE_LIB"
    JOBID=$(sbatch --parsable $SCRIPTS/STEP_07_annotate_breakpoints.sh \
        results/05_synteny_graph/breakpoints.bed \
        results/00_panSN_inputs/catfish_${TIER}.fa.gz \
        "$TE_ARG")
    echo "    Submitted job $JOBID"
    while squeue -j "$JOBID" &>/dev/null; do sleep 60; done
fi

# ---- STAGE 8: cross-species stats ----
run_stage 8 "Cross-species satellite stats" \
    "python3 $SCRIPTS/STEP_08_cross_species_stats.py \
        --annotation-dir results/07_bp_annotation/ \
        --manifest config/species_manifest.tsv \
        --out results/08_cross_species_stats/" \
    "results/08_cross_species_stats/satellite_catalog.tsv"

# ---- STAGE 9: dual evidence (SLURM) ----
if [[ $STAGE_START -le 9 && $STAGE_END -ge 9 ]] && [[ -n "$PROTEOME" ]] && ! stage_done "results/09_dual_evidence/breakpoint_confidence.tsv"; then
    echo ""
    echo "=== Stage 9: dual-evidence validation (SLURM) ==="
    JOBID=$(sbatch --parsable $SCRIPTS/STEP_09_dual_evidence_validate.sh \
        results/05_synteny_graph/breakpoints.bed \
        results/00_panSN_inputs/catfish_${TIER}.fa.gz \
        "$PROTEOME")
    while squeue -j "$JOBID" &>/dev/null; do sleep 60; done
fi

# ---- STAGE 9c: flank coherence ----
if [[ $STAGE_START -le 9 && $STAGE_END -ge 9 ]] && [[ -n "$PROTEOME" ]] && ! stage_done "results/09c_wg_miniprot/.done"; then
    echo ""
    echo "=== Stage 9c1: whole-genome miniprot (SLURM array) ==="
    N_SP=$(tail -n +2 config/species_manifest.tsv | awk '$1 !~ /^#/' | wc -l)
    JOBID=$(sbatch --parsable --array=0-$((N_SP-1)) \
        $SCRIPTS/STEP_09c1_miniprot_wholegenome.sh "$PROTEOME" "$TIER")
    while squeue -j "$JOBID" &>/dev/null; do sleep 60; done
    touch results/09c_wg_miniprot/.done
fi

run_stage 9 "Flank coherence scoring" \
    "python3 $SCRIPTS/STEP_09c_flank_coherence.py \
        --bp-bed results/05_synteny_graph/breakpoints.bed \
        --wg-miniprot-dir results/09c_wg_miniprot/ \
        --manifest config/species_manifest.tsv \
        --out results/09c_flank_coherence/ \
        --flank-kb 500 --min-families 3 --max-gap-genes 10" \
    "results/09c_flank_coherence/flank_coherence.tsv"

# ---- STAGE 10: enriched graph ----
run_stage 10 "Enrich graph for Cytoscape" \
    "python3 $SCRIPTS/STEP_10_enrich_graph_for_cytoscape.py \
        --graph-pickle results/05_synteny_graph/graph.pickle \
        --bp-annotation results/07_bp_annotation/breakpoint_annotation_summary.tsv \
        --flank-coherence results/09c_flank_coherence/flank_coherence.tsv \
        --dual-evidence results/09_dual_evidence/breakpoint_confidence.tsv \
        --out results/10_enriched_graph/" \
    "results/10_enriched_graph/graph_enriched.gexf"

# ---- STAGE 11: tree polarization ----
run_stage 11 "Tree-based polarization" \
    "python3 $SCRIPTS/STEP_11_polarize_with_tree.py \
        --tree config/species_tree.nwk \
        --flank-details results/09c_flank_coherence/flank_match_details.tsv \
        --bp-coherence results/09c_flank_coherence/breakpoint_coherence.tsv \
        --out results/11_polarized_breakpoints/" \
    "results/11_polarized_breakpoints/polarized_breakpoints.tsv"

echo ""
echo "=============================================================="
echo "Pipeline complete. Final outputs:"
echo "  results/10_enriched_graph/graph_enriched.gexf    (open in Cytoscape)"
echo "  results/11_polarized_breakpoints/polarized_breakpoints.tsv"
echo "=============================================================="
