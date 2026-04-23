#!/usr/bin/env Rscript
# ============================================================================
# STEP_04_visualize_synteny.R
#   Produce two kinds of figures from the wfmash scaffold PAF:
#     1. Pairwise dotplots (one per species pair) — quick visual check
#     2. Linear synteny ribbons (Kuang et al. Figure 1B style) — publication plot
#
# Usage:
#   Rscript STEP_04_visualize_synteny.R <tier>
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
tier <- if (length(args) >= 1) args[1] else "core"

paf_path <- sprintf("results/01_wfmash/catfish_%s.scaffolds.paf", tier)
bp_path  <- "results/02_breakpoints/breakpoints_all.tsv"
outdir   <- "results/04_figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# Load PAF
# ----------------------------------------------------------------------------
paf_cols <- c("qname","qlen","qstart","qend","strand",
              "tname","tlen","tstart","tend",
              "matches","aln_len","mapq")

paf <- read_tsv(paf_path, col_names = paf_cols,
                col_types = cols(.default = col_character(),
                                 qlen = col_integer(), qstart = col_integer(), qend = col_integer(),
                                 tlen = col_integer(), tstart = col_integer(), tend = col_integer(),
                                 matches = col_integer(), aln_len = col_integer(), mapq = col_integer()),
                col_select = 1:12)

# Split PanSN names
paf <- paf %>%
  separate(qname, into = c("qsp", "qhap", "qchrom"), sep = "#", remove = FALSE, extra = "merge") %>%
  separate(tname, into = c("tsp", "thap", "tchrom"), sep = "#", remove = FALSE, extra = "merge")

cat(sprintf("Loaded %d scaffold records\n", nrow(paf)))
cat(sprintf("Species pairs: %d\n",
            paf %>% distinct(qsp, tsp) %>% nrow()))

# ----------------------------------------------------------------------------
# Breakpoints for overlay
# ----------------------------------------------------------------------------
bp <- if (file.exists(bp_path)) {
  read_tsv(bp_path, show_col_types = FALSE)
} else {
  tibble()
}

# ----------------------------------------------------------------------------
# Chromosome ordering helper — numeric sort of e.g. LG01, LG02, ..., LG10
# ----------------------------------------------------------------------------
chrom_order <- function(x) {
  # extract trailing number if any; fall back to alphabetic
  nums <- as.integer(str_extract(x, "\\d+$"))
  x[order(nums, x, na.last = TRUE)]
}

# ============================================================================
# FIGURE 1: Pairwise dotplots
# ============================================================================
make_dotplot <- function(df, qsp_i, tsp_i) {
  d <- df %>% filter(qsp == qsp_i, tsp == tsp_i)
  if (nrow(d) == 0) return(NULL)

  q_levels <- chrom_order(unique(d$qchrom))
  t_levels <- chrom_order(unique(d$tchrom))
  d$qchrom <- factor(d$qchrom, levels = q_levels)
  d$tchrom <- factor(d$tchrom, levels = t_levels)

  ggplot(d) +
    geom_segment(aes(x = qstart, xend = qend,
                     y = tstart, yend = tend,
                     color = strand),
                 linewidth = 0.6) +
    facet_grid(tchrom ~ qchrom, scales = "free", space = "free", switch = "both") +
    scale_color_manual(values = c("+" = "#1f77b4", "-" = "#d62728")) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey80"),
          panel.spacing = unit(0.05, "lines"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.text.x = element_text(size = 6, angle = 90),
          strip.text.y.left = element_text(size = 6, angle = 0),
          strip.placement = "outside",
          legend.position = "top") +
    labs(title = sprintf("%s (query) vs %s (target)", qsp_i, tsp_i),
         x = qsp_i, y = tsp_i,
         color = "Strand")
}

pairs <- paf %>% distinct(qsp, tsp) %>% filter(qsp < tsp)  # dedupe symmetric pairs

for (i in seq_len(nrow(pairs))) {
  qsp_i <- pairs$qsp[i]; tsp_i <- pairs$tsp[i]
  p <- make_dotplot(paf, qsp_i, tsp_i)
  if (!is.null(p)) {
    outfile <- file.path(outdir, sprintf("dotplot_%s_vs_%s.pdf", qsp_i, tsp_i))
    ggsave(outfile, p, width = 12, height = 12)
    cat(sprintf("  Dotplot: %s\n", outfile))
  }
}

# ============================================================================
# FIGURE 2: Linear synteny ribbons (Kuang Figure 1B style)
# Stacks species horizontally; ribbons connect homologous regions.
# ============================================================================
make_ribbon_plot <- function(paf_df, species_order) {
  # Build chromosome layout — one row per species, chroms laid out side by side
  chrom_lens <- paf_df %>%
    pivot_longer(c(qsp, tsp), names_to = "role", values_to = "sp") %>%
    mutate(chrom = ifelse(role == "qsp", qchrom, tchrom),
           len   = ifelse(role == "qsp", qlen,   tlen)) %>%
    distinct(sp, chrom, len) %>%
    filter(!is.na(chrom)) %>%
    group_by(sp) %>%
    arrange(match(chrom, chrom_order(chrom))) %>%
    mutate(chrom_cumstart = cumsum(lag(len, default = 0)) + (row_number() - 1) * 5e6) %>%
    ungroup()

  layout <- chrom_lens %>%
    mutate(chrom_start_abs = chrom_cumstart,
           chrom_end_abs   = chrom_cumstart + len,
           y = match(sp, species_order))

  # Build ribbons: each PAF row becomes a polygon between (qsp row) and (tsp row)
  ribbons <- paf_df %>%
    inner_join(layout %>% select(qsp = sp, qchrom = chrom, q_offset = chrom_start_abs),
               by = c("qsp", "qchrom")) %>%
    inner_join(layout %>% select(tsp = sp, tchrom = chrom, t_offset = chrom_start_abs),
               by = c("tsp", "tchrom")) %>%
    mutate(q_abs_start = qstart + q_offset,
           q_abs_end   = qend   + q_offset,
           t_abs_start = tstart + t_offset,
           t_abs_end   = tend   + t_offset,
           qy = match(qsp, species_order),
           ty = match(tsp, species_order))

  # Plot
  p <- ggplot() +
    # Chromosome bars
    geom_rect(data = layout,
              aes(xmin = chrom_start_abs, xmax = chrom_end_abs,
                  ymin = y - 0.1, ymax = y + 0.1,
                  fill = chrom),
              color = "black", linewidth = 0.2) +
    # Ribbons (simplified: straight segments between block midpoints)
    geom_segment(data = ribbons,
                 aes(x = (q_abs_start + q_abs_end)/2, xend = (t_abs_start + t_abs_end)/2,
                     y = qy, yend = ty,
                     color = strand),
                 alpha = 0.15, linewidth = 0.3) +
    scale_fill_viridis_d(guide = "none") +
    scale_color_manual(values = c("+" = "grey30", "-" = "#d62728")) +
    scale_y_continuous(breaks = seq_along(species_order), labels = species_order) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(x = NULL, y = NULL, color = "Strand",
         title = "Genome-scale synteny across catfish")

  p
}

species_order <- sort(unique(c(paf$qsp, paf$tsp)))
if (length(species_order) >= 2) {
  p_ribbon <- make_ribbon_plot(paf, species_order)
  ggsave(file.path(outdir, "synteny_ribbons.pdf"),
         p_ribbon, width = 16, height = 2 + length(species_order))
  cat(sprintf("  Ribbon plot: %s/synteny_ribbons.pdf\n", outdir))
}

# ============================================================================
# FIGURE 3: Breakpoint summary barplot
# ============================================================================
if (nrow(bp) > 0) {
  bp_summary <- bp %>%
    count(query_species, target_species, event_type) %>%
    mutate(pair = paste(query_species, "vs", target_species))

  p_bp <- ggplot(bp_summary, aes(x = pair, y = n, fill = event_type)) +
    geom_col(position = "dodge") +
    coord_flip() +
    scale_fill_manual(values = c("FUSION_in_query"       = "#1f77b4",
                                 "INVERSION_breakpoint"  = "#ff7f0e",
                                 "NONCOLLINEAR"          = "grey60")) +
    theme_minimal() +
    labs(title = "Breakpoint events per species pair",
         x = NULL, y = "Number of events", fill = "Event type")

  ggsave(file.path(outdir, "breakpoint_summary.pdf"), p_bp,
         width = 10, height = 1 + 0.3 * nrow(bp_summary))
  cat(sprintf("  Breakpoint summary: %s/breakpoint_summary.pdf\n", outdir))
}

cat("[STEP_04] Done.\n")
