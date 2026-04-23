# MANUSCRIPT_SKELETON.md

Methods section scaffold for the catfish synteny/breakpoint manuscript.
Everything in `<FILL IN: ...>` refers to specific output files — replace
with the actual numbers/values after running the pipeline.

**Do not write prose interpretations before you have results.**
**Do not over-promise in Methods; describe what you did, not what you hoped for.**

---

## Methods

### Genome assemblies and data sources

We analyzed chromosome-level assemblies of `<FILL IN: N>` catfish species:
`<FILL IN: list species with genus + assembly accession from config/ncbi_accessions.tsv>`.
Assemblies were obtained from NCBI GenBank/RefSeq (accessed
`<FILL IN: date>`) using HTTPS downloads with md5 verification
(`download_ncbi_genome.sh` in the catfish-synteny-toolkit). The NCBI
datasets CLI was not used due to documented reliability issues with
chained batch downloads.

Per-species assembly statistics are summarized in
**Supplementary Table S1** (see `results/00_prepared_genomes/<sp>/<sp>.report.txt`
for each species).

### Genome preparation

Each assembly was processed through the ScaffoldKit Module 1 pipeline
(STEP_00b_prepare_genomes.sh): FASTA normalization with Picard
NormalizeFasta (line length 80 bp, headers truncated at whitespace);
scaffold size tabulation with seqkit fx2tab; gap detection with detgaps;
AGP generation and contig splitting with RagTag splitasm; telomere
detection with Quartet TeloExplorer; and removal of contigs shorter than
5 Mb using seqkit to retain chromosome-level scaffolds only. MD5
checksums of all processed files are provided in
`results/00_prepared_genomes/<sp>/<sp>.md5`.

### All-vs-all homology mapping (Stage A)

Processed genomes were concatenated using the PanSN naming convention
(`<species>#1#<chromosome>`) and indexed with samtools faidx. All-vs-all
homology mapping was performed with wfmash
(version `<FILL IN: cat versions/wfmash.txt>`), which combines MashMap3
minmer sketching with WFA base-level alignment. Key parameters:

- **Segment length (-s)**: `<FILL IN: 50000 for core/clarias tiers,
  30000 for all tier, 20000 for deep tier>`
- **Minimum block length (-l)**: 3 × segment length
- **Mapping identity (-p)**: `<FILL IN: ani50-3 for core tier (~90% floor),
  ani50-5 for clarias, ani50-15 for all, -p 70 for deep outgroup>`
- **Scaffold mass (-S)**: 100 kb (or `<FILL IN>` for tier)
- **Frequent k-mer filter (-F)**: 0.00005 (aggressive, appropriate for
  repeat-rich fish genomes)
- **PanSN skip (-Y '#')**: intra-species mappings excluded
- **Filter (-o, -n 1)**: one-to-one plane-sweep filter

Intra-species mappings were excluded via the `-Y '#'` PanSN delimiter.
Intermediate scaffold PAF (`--scaffold-out`) was retained as the primary
input to downstream breakpoint calling.

### Synteny graph construction (Stage B)

From the scaffold PAF, we constructed an undirected multi-graph
(STEP_05_build_synteny_graph.py) where nodes are genomic segments
bounded by inferred breakpoints and edges are wfmash homology blocks
weighted by block length and sequence identity. Breakpoints were
identified as positions where adjacent blocks on the same query
chromosome mapped to different target chromosomes (inferred fusion/
fission signals) or switched strand (inferred inversion breakpoints).
Breakpoints within 500 kb of each other were merged. Blocks shorter
than 200 kb were excluded to minimize transposon-driven noise.

The resulting graph contained `<FILL IN: N_nodes from
results/05_synteny_graph/node_summary.tsv>` nodes and
`<FILL IN: N_edges>` homology edges connecting
`<FILL IN: N_species>` species.

### Breakpoint annotation (Stage C)

Each breakpoint window (±`<FILL IN: 100>` kb) was annotated using a
panel of established tools (STEP_07_annotate_breakpoints.sh):

- **Tandem repeats and satellites**: TRF version
  `<FILL IN: cat versions/trf.txt>` with recommended parameters
  2 7 7 80 10 50 500. Repeats with period ≥ 100 bp were classified as
  satellite-class; period 1–6 bp as microsatellites.
- **Inverted repeats and dyad symmetry**: IRF version
  `<FILL IN: cat versions/irf.txt>` with parameters 2 3 5 80 10 40
  500000 10000.
- **Segmental duplications**: minimap2 version
  `<FILL IN: cat versions/minimap2.txt>` in self-alignment mode
  (-X -cx asm20); non-identity hits ≥ 1 kb were retained.
- **Transposable elements**: RepeatMasker with a catfish-specific
  library (`<FILL IN: describe TE library source — e.g. from the LG28
  inversion project, built with EDTA or custom>`).

### Cross-species satellite analysis (Stage D)

Satellite consensus sequences identified by TRF at breakpoints were
queried against each species' full genome using minimap2 (map-ont
preset). Hits were tabulated per chromosome per species to generate an
enrichment matrix analogous to the analysis in Kuang et al. 2026
(STEP_08_cross_species_stats.py). A satellite was considered enriched
at a breakpoint chromosome if its hit count there exceeded the
genome-wide per-chromosome mean (chi-square test).

### Multi-evidence breakpoint validation (Stage E)

Three independent evidence types were integrated per breakpoint:

1. **k-mer homology robustness** (STEP_09_dual_evidence_validate.sh):
   mashmap was run in a multi-pi sweep (pi = 70, 80, 85, 90, 95) on
   breakpoint windows; each 5 kb tile was scored by the highest pi
   threshold at which it retained homology to another window.

2. **Per-protein support** (STEP_09b_score_dual_evidence.py): miniprot
   version `<FILL IN: cat versions/miniprot.txt>` with `<FILL IN:
   reference proteome: e.g. C. gariepinus predicted proteins or
   zebrafish RefSeq>` aligned to breakpoint windows; each tile scored
   by the presence of ≥1 miniprot mRNA feature.

3. **Protein-family order preservation** (STEP_09c_flank_coherence.py):
   miniprot was run genome-wide on each species with lenient
   parameters (`--outs=0.5`). For each breakpoint flank (±500 kb), the
   ordered list of matched protein families was extracted. A flank
   was scored as showing STRONG_ORTHOLOGY if ≥60% of families appeared
   in the same order (allowing up to 10 intervening non-matching
   genes) in at least one other species' genome — a collect-then-refine
   adaptation of the hamburger clustering algorithm (djw533).

Breakpoints were classified as HIGH confidence when both flanks showed
STRONG or PARTIAL orthology; NOVELTY when exactly one flank was
orthologous (the other being lineage-specific or gene-poor); and LOW
or UNINFORMATIVE otherwise.

### Breakpoint polarization using species tree (Stage F)

A BUSCO single-copy supermatrix phylogeny
(`<FILL IN: method — e.g. IQ-TREE with LG+G, N single-copy orthologs>`)
was used as the reference species tree (config/species_tree.nwk). For
each breakpoint, the set of species sharing each flank's gene order was
mapped onto the tree (STEP_11_polarize_with_tree.py). A breakpoint was
polarized by parsimony: if one flank's sharer set was a strict subset
of the other's, the rearrangement was inferred to have occurred on the
branch leading to the species carrying the derived state. Ambiguous
configurations (non-nested sharer sets) were reported as such and not
forced to a direction.

### Reproducibility

All scripts, configuration files, and the exact command-line invocations
are archived in the catfish-synteny-toolkit (version `<FILL IN: git
commit hash or release tag>`). Tool versions were captured at run-time
(see Supplementary Methods). Input genome MD5 checksums are in
`results/00_prepared_genomes/*/md5`. The enriched synteny graph
(`results/10_enriched_graph/graph_enriched.gexf`) is provided for
interactive exploration in Cytoscape or Gephi.

---

## Results section — sections to write, not sections to draft

**Do not draft these until you have data.** Use as section headers:

1. Genome assembly characteristics
   - Table of assembly stats (from `results/00_prepared_genomes/*/report.txt`)
   - Figure: cumulative assembly size, gap distributions

2. Genome-wide synteny between catfish species
   - Figure: synteny ribbons (`results/05_synteny_graph/figures/synteny_ribbons.pdf`)
   - Figure: pairwise synteny Mb matrix (`pair_heatmap.pdf`)

3. Structural rearrangements between Clarias species
   - Count and categorize breakpoints (`breakpoint_coherence.tsv`)
   - Karyotype-scale summary: which chromosomes are rearranged

4. Mechanistic signatures at breakpoints
   - IRF palindrome enrichment (`breakpoint_annotation_summary.tsv`)
   - Satellite enrichment (`Fig2_equivalent.pdf`)
   - TE content comparison vs genome background

5. Multi-evidence validation of breakpoints
   - Distribution of confidence classes
   - NOVELTY cases: describe the 1–3 most striking examples (but ONLY
     if the pipeline actually produced them)

6. Tree-polarized rearrangement history
   - Breakpoints polarized by T. rosablanca outgroup
   - (Honest limitation: polarization depth limited by outgroup
     divergence and absence of divergence-time calibration)

---

## What to NOT write in Methods or Results

- Do not claim ancestral karyotype reconstruction unless you actually
  ran it and the tree topology supported unique polarization.
- Do not estimate divergence times — you don't have fossils or a
  molecular clock calibration. Say so explicitly.
- Do not describe Fisher's exact tests, permutations, or other
  statistics unless the scripts actually compute them. Current Stage D
  uses a chi-square approximation; don't upgrade the description.
- Do not claim biological mechanism (e.g., "palindrome-mediated
  rearrangement") unless you see a clear enrichment signal. If absent,
  report absence honestly.

---

## Reviewer defense notes (for later)

When you receive reviewer comments, the likely asks will be:

- "Why not WGDI / MCScanX?" → See DEPENDENCIES.md section on excluded
  tools. wfmash is more modern and better suited to vertebrate pairwise
  synteny.
- "Why not Hi-C / 3D architecture?" → Out of scope; data not available.
  Cite as limitation.
- "How do you distinguish true rearrangement from assembly artifacts?"
  → Multi-evidence validation (Stage E) with three independent layers.
- "Sample size is small (N species)." → Acknowledge; note that
  chromosome-level catfish genomes are rare and this is the largest
  comparative set currently possible.
- "Polarization assumes tree is correct." → Acknowledged; BUSCO
  supermatrix tree is cited. Alternative topologies tested if you want
  to defend further.
