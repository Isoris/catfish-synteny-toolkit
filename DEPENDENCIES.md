# DEPENDENCIES.md

Versioned tool matrix for catfish-synteny-toolkit. Pin these versions when
publishing. BibTeX entries at the bottom.

## Core alignment tools

| Tool | Version pin | Purpose | Citation key |
|------|-------------|---------|--------------|
| wfmash | ≥ 0.14.0 | All-vs-all homology mapping (Stage A) | @guarracino2024 |
| MashMap3 | bundled with wfmash ≥ 0.14 | k-mer-based approximate mapping | @kille2023 |
| minimap2 | ≥ 2.24 | Self-alignment for SDs; local refinement | @li2018minimap2 |
| miniprot | ≥ 0.12 | Protein-to-genome alignment | @li2023miniprot |
| samtools | ≥ 1.17 | FASTA indexing | @danecek2021 |
| bcftools | ≥ 1.17 | (optional, for downstream) | @danecek2021 |
| htslib | ≥ 1.17 | bgzip compression | @danecek2021 |

## Repeat / structural annotation

| Tool | Version pin | Purpose | Citation key |
|------|-------------|---------|--------------|
| TRF | ≥ 4.09 | Tandem repeats, satellites, microsatellites | @benson1999trf |
| IRF | ≥ 3.08 | Inverted repeats (dyad symmetry) | @warburton2004irf |
| RepeatMasker | ≥ 4.1.5 (optional) | TE annotation with user library | @smit2013repeatmasker |
| detgaps | (part of ScaffoldKit) | Gap detection | N/A |
| seqkit | ≥ 2.5.0 | FASTA manipulation + size filter | @shen2016seqkit |

## Genome preparation

| Tool | Version pin | Purpose | Citation key |
|------|-------------|---------|--------------|
| Picard | ≥ 3.0 | NormalizeFasta | @broadinstitute_picard |
| RagTag | ≥ 2.1.0 | splitasm → contig AGP | @alonge2022ragtag |
| Quartet TeloExplorer | ≥ 1.1.4 | Telomere detection | @lin2023quartet |
| ScaffoldKit Module 1 | custom (project-internal) | Per-species prep orchestration | N/A |

## Graph / visualization / analysis

| Tool | Version pin | Purpose | Citation key |
|------|-------------|---------|--------------|
| networkx | ≥ 3.0 | Synteny graph construction | @hagberg2008networkx |
| matplotlib | ≥ 3.7 | All Python figures | @hunter2007matplotlib |
| numpy | ≥ 1.24 | Array math | @harris2020numpy |
| ete3 | ≥ 3.1.3 | Species tree parsing for polarization | @huerta2016ete3 |
| Cytoscape | ≥ 3.10 | Interactive graph visualization (user-side) | @shannon2003cytoscape |
| Gephi | ≥ 0.10 (alternative to Cytoscape) | GEXF visualization | @bastian2009gephi |
| R | ≥ 4.3 | Optional R dotplots (STEP_04) | N/A |
| ggplot2, dplyr, tidyr | CRAN | R visualization | @wickham2016ggplot2 |

## Language runtimes

| Runtime | Version pin | Used in |
|---------|-------------|---------|
| Python | ≥ 3.10 | STEP_02, 05, 06, 08, 09b, 09c, 10, 11 |
| Bash | ≥ 4.4 | All `.sh` scripts |
| R | ≥ 4.3 | STEP_04 (optional) |

## Data sources

| Source | URL | Accession scheme | Access method |
|--------|-----|------------------|---------------|
| NCBI GenBank/RefSeq | https://ftp.ncbi.nlm.nih.gov/genomes/all/ | GCA_/GCF_ | HTTPS wget (toolkit's download_ncbi_genome.sh) |
| NCBI datasets API | https://www.ncbi.nlm.nih.gov/datasets/ | GCA_/GCF_ | Not used (unreliable for batch) |

**Note**: NCBI rsync support retires 2026-06-01. This toolkit does not depend
on rsync.

## Capturing versions at run-time

Before a publication-ready run, capture the exact versions used:

```bash
# Add to the top of your SLURM job or a pre-run script
mkdir -p versions
wfmash --version > versions/wfmash.txt
mashmap --version > versions/mashmap.txt 2>&1
minimap2 --version > versions/minimap2.txt
miniprot --version > versions/miniprot.txt
samtools --version | head -1 > versions/samtools.txt
bcftools --version | head -1 > versions/bcftools.txt
trf 2>&1 | head -2 > versions/trf.txt || true
irf 2>&1 | head -2 > versions/irf.txt || true
seqkit version > versions/seqkit.txt
RepeatMasker -v 2>/dev/null | head -1 > versions/repeatmasker.txt || true
python3 --version > versions/python.txt
python3 -c "import networkx, matplotlib, numpy, ete3; \
    print(f'networkx={networkx.__version__}'); \
    print(f'matplotlib={matplotlib.__version__}'); \
    print(f'numpy={numpy.__version__}'); \
    print(f'ete3={ete3.__version__}')" > versions/python_packages.txt
```

Include `versions/` in your supplementary materials. This gives reviewers
byte-exact reproducibility without depending on conda-env YAMLs that
drift over time.

## BibTeX entries

```bibtex
@article{guarracino2024,
  title   = {Building pangenome graphs},
  author  = {Guarracino, Andrea and Heumos, Simon and Nahnsen, Sven and
             Prins, Pjotr and Garrison, Erik},
  journal = {Bioinformatics},
  year    = {2024},
  note    = {wfmash component},
}

@article{kille2023,
  title   = {Minmers are a generalization of minimizers that enable unbiased
             local Jaccard estimation},
  author  = {Kille, Bryce and Garrison, Erik and Treangen, Todd J and
             Phillippy, Adam M},
  journal = {Bioinformatics},
  year    = {2023},
  volume  = {39},
  number  = {9},
  pages   = {btad512},
  doi     = {10.1093/bioinformatics/btad512},
}

@article{li2018minimap2,
  title   = {Minimap2: pairwise alignment for nucleotide sequences},
  author  = {Li, Heng},
  journal = {Bioinformatics},
  year    = {2018},
  volume  = {34},
  number  = {18},
  pages   = {3094--3100},
  doi     = {10.1093/bioinformatics/bty191},
}

@article{li2023miniprot,
  title   = {Protein-to-genome alignment with miniprot},
  author  = {Li, Heng},
  journal = {Bioinformatics},
  year    = {2023},
  volume  = {39},
  number  = {1},
  pages   = {btad014},
  doi     = {10.1093/bioinformatics/btad014},
}

@article{danecek2021,
  title   = {Twelve years of {SAMtools} and {BCFtools}},
  author  = {Danecek, Petr and Bonfield, James K and Liddle, Jennifer and
             Marshall, John and Ohan, Valeriu and Pollard, Martin O and
             Whitwham, Andrew and Keane, Thomas and McCarthy, Shane A and
             Davies, Robert M and Li, Heng},
  journal = {GigaScience},
  year    = {2021},
  volume  = {10},
  number  = {2},
  pages   = {giab008},
  doi     = {10.1093/gigascience/giab008},
}

@article{benson1999trf,
  title   = {Tandem repeats finder: a program to analyze {DNA} sequences},
  author  = {Benson, Gary},
  journal = {Nucleic Acids Research},
  year    = {1999},
  volume  = {27},
  number  = {2},
  pages   = {573--580},
  doi     = {10.1093/nar/27.2.573},
}

@article{warburton2004irf,
  title   = {Inverted repeat structure of the human genome: the
             X-chromosome contains a preponderance of large, highly
             homologous inverted repeats},
  author  = {Warburton, Peter E and Giordano, J and Cheung, F and
             Gelfand, Y and Benson, Gary},
  journal = {Genome Research},
  year    = {2004},
  volume  = {14},
  number  = {10a},
  pages   = {1861--1869},
  doi     = {10.1101/gr.2542904},
}

@misc{smit2013repeatmasker,
  title        = {{RepeatMasker Open-4.0}},
  author       = {Smit, AFA and Hubley, R and Green, P},
  year         = {2013--2015},
  howpublished = {\url{http://www.repeatmasker.org}},
}

@article{shen2016seqkit,
  title   = {{SeqKit}: A cross-platform and ultrafast toolkit for
             {FASTA/Q} file manipulation},
  author  = {Shen, Wei and Le, Shuai and Li, Yan and Hu, Fuquan},
  journal = {PLOS ONE},
  year    = {2016},
  volume  = {11},
  number  = {10},
  pages   = {e0163962},
  doi     = {10.1371/journal.pone.0163962},
}

@misc{broadinstitute_picard,
  title        = {Picard Tools},
  author       = {{Broad Institute}},
  howpublished = {\url{http://broadinstitute.github.io/picard/}},
  year         = {2019},
}

@article{alonge2022ragtag,
  title   = {Automated assembly scaffolding using {RagTag} elevates a new
             tomato system for high-throughput genome editing},
  author  = {Alonge, Michael and Lebeigle, Ludivine and Kirsche, Melanie
             and Jenike, Katie and Ou, Shujun and Aganezov, Sergey and
             Wang, Xingang and Lippman, Zachary B and Schatz, Michael C
             and Soyk, Sebastian},
  journal = {Genome Biology},
  year    = {2022},
  volume  = {23},
  pages   = {258},
  doi     = {10.1186/s13059-022-02823-7},
}

@article{lin2023quartet,
  title   = {quarTeT: a telomere-to-telomere toolkit for gap-free genome
             assembly and centromeric repeat identification},
  author  = {Lin, Yun and Ye, Changlin and Li, Xinru and Chen, Qiang and
             Wu, Yinyue and Zhang, Fang and Pan, Runze and Zhang, Saining
             and Chen, Shuaibing and Wang, Shuo and Lv, Shuo and Pan, Yang
             and Huang, Ying and Lu, Yangshuo and Cao, Shuqi and Wang,
             Yi},
  journal = {Horticulture Research},
  year    = {2023},
  volume  = {10},
  pages   = {uhad127},
  doi     = {10.1093/hr/uhad127},
}

@inproceedings{hagberg2008networkx,
  title     = {Exploring network structure, dynamics, and function using
               {NetworkX}},
  author    = {Hagberg, Aric A and Schult, Daniel A and Swart, Pieter J},
  booktitle = {Proceedings of the 7th Python in Science Conference},
  year      = {2008},
  pages     = {11--15},
}

@article{hunter2007matplotlib,
  title   = {Matplotlib: a {2D} graphics environment},
  author  = {Hunter, John D},
  journal = {Computing in Science \& Engineering},
  year    = {2007},
  volume  = {9},
  number  = {3},
  pages   = {90--95},
  doi     = {10.1109/MCSE.2007.55},
}

@article{harris2020numpy,
  title   = {Array programming with {NumPy}},
  author  = {Harris, Charles R and Millman, K Jarrod and van der Walt,
             Stéfan J and Gommers, Ralf and Virtanen, Pauli and
             Cournapeau, David and others},
  journal = {Nature},
  year    = {2020},
  volume  = {585},
  number  = {7825},
  pages   = {357--362},
  doi     = {10.1038/s41586-020-2649-2},
}

@article{huerta2016ete3,
  title   = {{ETE 3}: reconstruction, analysis, and visualization of
             phylogenomic data},
  author  = {Huerta-Cepas, Jaime and Serra, François and Bork, Peer},
  journal = {Molecular Biology and Evolution},
  year    = {2016},
  volume  = {33},
  number  = {6},
  pages   = {1635--1638},
  doi     = {10.1093/molbev/msw046},
}

@article{shannon2003cytoscape,
  title   = {Cytoscape: a software environment for integrated models of
             biomolecular interaction networks},
  author  = {Shannon, Paul and Markiel, Andrew and Ozier, Owen and
             Baliga, Nitin S and Wang, Jonathan T and Ramage, Daniel and
             Amin, Nada and Schwikowski, Benno and Ideker, Trey},
  journal = {Genome Research},
  year    = {2003},
  volume  = {13},
  number  = {11},
  pages   = {2498--2504},
  doi     = {10.1101/gr.1239303},
}

@inproceedings{bastian2009gephi,
  title     = {Gephi: an open source software for exploring and manipulating
               networks},
  author    = {Bastian, Mathieu and Heymann, Sebastien and Jacomy, Mathieu},
  booktitle = {Proceedings of the International AAAI Conference on Web and
               Social Media},
  year      = {2009},
  volume    = {3},
  number    = {1},
  pages     = {361--362},
}

@book{wickham2016ggplot2,
  title     = {ggplot2: Elegant Graphics for Data Analysis},
  author    = {Wickham, Hadley},
  publisher = {Springer},
  year      = {2016},
  edition   = {2nd},
}
```

## Adjacent tools referenced but not used (for completeness)

Tools that appeared in discussions but were **deliberately excluded**:

- **WGDI** / **MCScanX** — plant-WGD-oriented collinearity tools, slower and
  weaker than wfmash for vertebrate pairwise synteny. Not used.
- **TideCluster** / **TRASH** — satellite identification tools with limited
  validation; TRF + IRF chosen instead for ~25 years of community use.
- **SyRI** — structural variant identifier; overlaps wfmash+custom-parsing
  without clear advantage for this project.
- **Hamburger** (djw533) — inspired the flank coherence algorithm in STEP_09c
  but the bacterial-specific code was not reused directly.
