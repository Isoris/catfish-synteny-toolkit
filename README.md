# catfish-synteny-toolkit

Graph-based structural rearrangement detection and validation across catfish
genomes. Uses only established tools. No MCScanX/WGDI.

## One-command run

```bash
conda activate assembly
bash run_pipeline.sh core --stage-end 6        # fast scoping, ~2h
# look at results/05_synteny_graph/figures/  -- decide if signal is real
bash run_pipeline.sh clarias \
    --stage-start 7 \
    --proteome /path/to/cgar.faa \
    --te-lib /path/to/catfish_TE.fa           # full pipeline, ~12-24h
```

## Tools used (all established)

| Task | Tool | Citation |
|------|------|----------|
| All-vs-all homology | wfmash (MashMap3 + WFA) | Garrison et al. 2023 |
| Protein-to-genome | miniprot | Li 2023 Bioinformatics |
| Tandem repeats | TRF | Benson 1999 NAR |
| Palindromes / dyad | IRF | Warburton et al. 2004 |
| Self-alignments | minimap2 | Li 2018 Bioinformatics |
| TE annotation | RepeatMasker + user lib | Smit et al. |
| Graph | networkx + custom | Hagberg et al. 2008 |
| Tree parsing | ete3 | Huerta-Cepas et al. 2016 |

Deliberately NOT used: TideCluster, TRASH, SyRI, MCScanX, WGDI, GRIMM, MGR.

## Pipeline structure

```
Stage A (wfmash) → PAF → Stage B (graph) → breakpoints.bed
                                              ↓
                             ┌────────────────┼─────────────────┐
                             ↓                ↓                 ↓
                     Stage C (annotate) Stage D (stats)  Stage E (validate)
                     TRF, IRF, SDs, TEs  satellite matrix  3 evidence layers
                             └────────────────┼─────────────────┘
                                              ↓
                                       Stage F (integrate)
                                       enriched GEXF + tree polarization
```

## Species manifest tiers

| Tier | Use for | ANI strategy |
|------|---------|-------------|
| `core` | Cgar + Cmac only | ani50-3 (~90% floor) |
| `clarias` | core + other Clarias | ani50-5 |
| `all` | + Silurus/Tachysurus | ani50-15 |
| `deep` | + T. rosablanca | -p 70 (100 My divergence) |

Run tiers sequentially, not mixed in one pass.

## Key parameters

- **wfmash `-F 0.00005`**: aggressive frequent-kmer filter for repeat-rich fish
- **wfmash `-S 100000`**: scaffold mass — keep synteny blocks ≥100 kb
- **STEP_09c `--flank-kb 500`**: 500 kb flanks (100 kb too small for gene-poor regions)
- **STEP_09c `--min-families 3`**: hamburger-style clustering threshold

## Outputs that matter

| File | Contains |
|------|----------|
| `results/05_synteny_graph/figures/synteny_ribbons.pdf` | Kuang Fig 1B style |
| `results/05_synteny_graph/breakpoints.bed` | Input to all downstream |
| `results/07_bp_annotation/breakpoint_annotation_summary.tsv` | TRF/IRF/TE/SD per bp |
| `results/08_cross_species_stats/Fig2_equivalent.pdf` | Kuang Fig 2D-F style |
| `results/09c_flank_coherence/breakpoint_coherence.tsv` | Protein-order evidence |
| `results/10_enriched_graph/graph_enriched.gexf` | Cytoscape-ready with all scores |
| `results/11_polarized_breakpoints/polarized_breakpoints.tsv` | Tree-polarized direction |

## Decision checkpoints

**After Stage A/B (2 hours):** if no clean breakpoint signal between Cgar and
Cmac in the ribbon plot, STOP. Refocus on the inversion paper.

**After Stage C (full pipeline):** if breakpoints lack IRF palindromes and lack
satellite enrichment, you have a "we found rearrangements" paper, not a
"we found a mechanism" paper. Decide scope accordingly.

## Documented limitations

- wfmash scaffold filter misses sub-100 kb rearrangements
- Gene-poor flanks (desert regions) get TOO_FEW_GENES in STEP_09c — try
  `--flank-kb 500` or `1000`
- T. rosablanca at 100 My: protein anchors only, no k-mer synteny
- TE annotation quality depends entirely on user-provided library
- STEP_11 polarization assumes the BUSCO species tree is correct (user
  acknowledged the tree is "ugly but correct")

## Dependencies

```bash
conda install -c bioconda wfmash minimap2 miniprot samtools htslib trf
pip install --break-system-packages networkx matplotlib numpy ete3

# IRF not on bioconda — download from https://tandem.bu.edu/irf/
# RepeatMasker optional — only needed if using --te-lib
```

## See also

- `HANDOFF.md` — session audit and next-chat handoff
- `run_pipeline.sh --help` — orchestrator options
