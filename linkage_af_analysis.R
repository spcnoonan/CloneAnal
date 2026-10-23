#!/usr/bin/env Rscript
# ============================================================
# Eimeria cross: 50kb-binned LD vs allele-frequency-distortion overlay
#
# Inputs (from the PLINK pipeline):
#   eimeria_cross_ld.ld     -- pairwise r2 from `plink --r2`
#   eimeria_cross_freq.frq  -- per-site allele freq from `plink --freq`
#
# Output:
#   ld_freq_overlay.png     -- per-chromosome plot: mean r2 (50kb bins)
#                              vs allele-frequency skew, stacked
#   ld_freq_binned.csv      -- the underlying binned table, for ranking
#                              candidate blocks afterwards
# ============================================================

library(data.table)
library(ggplot2)
library(tidyverse)

# ------------------------------------------------------------
# STEP 0: User-set parameters
# ------------------------------------------------------------
ld_file   <- "vcftools_test.ld.reformat"
freq_file <- "clone_vars_freq.frq"
bin_size  <- 20000   # 50kb blocks, as requested

# ------------------------------------------------------------
# STEP 1: Load the LD table
# ------------------------------------------------------------
# plink --r2 output columns: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
# Each row is one SNP pair. To bin this into 20kb windows, we assign
# EACH SNP IN THE PAIR to its own genomic bin, then average. A pair
# spanning two different bins contributes to both bins' averages --
# this is intentional: it lets a strong long-range LD signal show up
# in every bin it touches, rather than disappearing because the two
# partner SNPs sit in different windows.
ld <- fread(ld_file)

# Defensive column name handling: plink sometimes pads column names
# with spaces depending on version; fread usually handles this, but
# we normalise just in case.
setnames(ld, trimws(names(ld)))

# ------------------------------------------------------------
# STEP 2: Load the allele frequency table
# ------------------------------------------------------------
# plink --freq output columns: CHR SNP A1 A2 MAF NCHROBS
# Here "MAF" is just the frequency of the A1 allele across the 18
# haploid clones -- since every site is biallelic between the two
# parents, this IS the parent-A vs parent-B allele ratio, not a
# population minor-allele frequency in the usual sense.
freq <- fread(freq_file)

setnames(freq, trimws(names(freq)))

# We need BP (position) for frequency sites, but plink --freq does not
# report position by default. Pull positions from the .bim file instead
# and merge on SNP id.
bim <- fread("all_vars.bim",
             col.names = c("CHR", "SNP", "CM", "BP", "A1_bim", "A2_bim"))

BP_loc = select(bim, BP)
freq = cbind(freq, BP_loc)

# ------------------------------------------------------------
# STEP 3: Assign genomic bins (20kb blocks) to both tables
# ------------------------------------------------------------
# Bin index = floor(position / bin_size); bin start coordinate used
# for plotting/joining is bin_index * bin_size.
freq[, bin := (BP %/% bin_size) * bin_size]

# For the LD table, assign each SNP of the pair to its own bin, then
# stack the two so both SNP_A's bin and SNP_B's bin get a contribution
# from this pair's r2. This is the "every pair counts toward every
# bin it touches" logic described in step 1.
ld_long <- rbindlist(list(
  ld[, .(CHR = CHR_A, BP = BP_A, R2)],
  ld[, .(CHR = CHR_B, BP = BP_B, R2)]
))
ld_long[, bin := (BP %/% bin_size) * bin_size]

# ------------------------------------------------------------
# STEP 4: Compute mean LD per 20kb bin, per chromosome
# ------------------------------------------------------------
ld_binned <- ld_long[, .(
  mean_r2 = mean(R2, na.rm = TRUE),
  n_pairs = .N
), by = .(CHR, bin)]

# ------------------------------------------------------------
# STEP 5: Compute mean allele frequency per 20kb bin, per chromosome
# ------------------------------------------------------------
freq_binned <- freq[, .(
  mean_af = mean(MAF, na.rm = TRUE),
  n_sites         = .N
), by = .(CHR, bin)]

# ------------------------------------------------------------
# STEP 6: Merge LD and distortion into one per-bin table
# ------------------------------------------------------------
# full outer join on CHR + bin, so bins with LD data but no freq data
# (or vice versa) are still retained and visible rather than silently
# dropped -- useful for spotting gaps in marker coverage.
binned <- merge(ld_binned, freq_binned, by = c("CHR", "bin"), all = TRUE)
setorder(binned, CHR, bin)

# A simple combined candidate score: bins that are high in BOTH LD and
# distortion rank highest. Each component is rescaled 0-1 within-genome
# first so that LD (0-1 already) and distortion (0-0.5) contribute
# comparably, then multiplied. This is a starting point, not a
# statistical test -- treat it as a way to rank/sort bins for visual
# follow-up, not a p-value.
binned[, mean_r2_norm        := mean_r2 / max(mean_r2, na.rm = TRUE)]
binned[, mean_af_norm:= mean_af / max(mean_af, na.rm = TRUE)]
binned[, candidate_score     := mean_r2_norm * mean_af_norm]

fwrite(binned, "ld_freq_binned.csv")

# ------------------------------------------------------------
# STEP 7: Plot -- LD and AF stacked per chromosome
# ------------------------------------------------------------
# Reshape to long format so LD and distortion can share an x-axis
# (genomic position) but sit in separate facet rows, faceted again
# by chromosome below that. This gives a Manhattan-style view where
# you visually scan for bins where both rows light up together.
plot_dt <- melt(
  binned,
  id.vars = c("CHR", "bin"),
  measure.vars = c("mean_r2", "mean_af"),
  variable.name = "metric",
  value.name = "value"
)

plot_dt[, metric := factor(metric,
                           levels = c("mean_r2", "mean_af"),
                           labels = c("Mean LD (r2)", "Mean allele-freq"))]

chr5_chr12 = filter(plot_dt, CHR == "HG994965.1" | CHR == "HG994972.1") %>%
  mutate(across('metric', str_replace, 'mean_r2', 'Mean LD (r2)')) %>%
  mutate(across('metric', str_replace, 'mean_af', 'Mean Allele Frequency')) %>%
  mutate(across('CHR', str_replace, 'HG994965.1', 'Chr5')) %>%
  mutate(across('CHR', str_replace, 'HG994972.1', 'Chr12'))

highlight <- tribble(
  ~CHR, ~metric, ~xmin, ~xmax,
  "Chr5", "Mean LD (r2)",          0.7, 1,
  "Chr5", "Mean Allele Frequency",       0.7, 1,
  "Chr12", "Mean LD (r2)",          1.5, 2,
  "Chr12", "Mean Allele Frequency",       1.5, 2
)

ggplot(chr5_chr12, aes(x = bin / 1e6, y = value)) +
  geom_line(linewidth = 0.4, colour = "grey40") +
  geom_point(size = 0.6, colour = "steelblue") +
  geom_rect(
    data = highlight,
    aes(xmin = xmin, xmax = xmax,
        ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "red",
    alpha = 0.2
  ) +
  facet_grid(metric ~ CHR, scales = "free_x", space = "free_x") +
  labs(
    x = "Position (Mb)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 15),
    panel.spacing.x = unit(2, "lines"),
    strip.background = element_rect(color="black", size=1, linetype="solid"),
    strip.text = element_text(size=15),
    axis.title.x = element_text(size=13)
    )

# For all chromosomes (one per row) - LD
exnuc = c("HG994976.1", "HG994977.1")
old <- c("HG994961.1", "HG994962.1", "HG994963.1","HG994964.1", "HG994965.1", "HG994966.1","HG994967.1", "HG994968.1", "HG994969.1","HG994970.1", "HG994971.1", "HG994972.1","HG994973.1", "HG994974.1", "HG994975.1")
new <- c("Chr1", "Chr2", "Chr3","Chr4", "Chr5", "Chr6","Chr7", "Chr8", "Chr9","Chr10", "Chr11", "Chr12","Chr13", "Chr14", "Chr15")
replacements <- c(
  "HG994961.1" = "Chr1",
  "HG994962.1" = "Chr2",
  "HG994963.1" = "Chr3",
  "HG994964.1" = "Chr4",
  "HG994965.1" = "Chr5",
  "HG994966.1" = "Chr6",
  "HG994967.1" = "Chr7",
  "HG994968.1" = "Chr8",
  "HG994969.1" = "Chr9",
  "HG994970.1" = "Chr10",
  "HG994971.1" = "Chr11",
  "HG994972.1" = "Chr12",
  "HG994973.1" = "Chr13",
  "HG994974.1" = "Chr14",
  "HG994975.1" = "Chr15"
)
plot_ld = filter(plot_dt, metric == "mean_r2") %>%
  filter(!(CHR %in% exnuc)) %>%
  mutate(across('metric', str_replace, 'mean_r2', 'Mean LD (r2)')) %>%
  mutate(across('metric', str_replace, 'mean_af', 'Mean Allele Frequency')) %>%
  mutate(across(CHR, ~ str_replace_all(.x, replacements))) %>%
  mutate(CHR = factor(CHR, levels=new))
  

ggplot(plot_ld, aes(x = bin / 1e6, y = value)) +
  geom_line(linewidth = 0.4, colour = "grey40") +
  geom_point(size = 0.6, colour = "steelblue") +
  facet_grid(rows = vars(CHR)) +
  labs(
    x = "Position (Mb)",
    y = "Mean LD (r2)\n"
  ) +
  theme_bw() +
  theme(
    strip.text.7 = element_text(size = 10),
    panel.spacing.y = unit(0.5, "lines"),
    strip.background = element_rect(color="black", size=1, linetype="solid"),
    strip.text = element_text(size=15),
    axis.title.x = element_text(size=13)
  )

# For all chromosomes (one per row) - LD
plot_af = filter(plot_dt, metric == "mean_af") %>%
  filter(!(CHR %in% exnuc)) %>%
  mutate(across('metric', str_replace, 'mean_r2', 'Mean LD (r2)')) %>%
  mutate(across('metric', str_replace, 'mean_af', 'Mean Allele Frequency')) %>%
  mutate(across(CHR, ~ str_replace_all(.x, replacements))) %>%
  mutate(CHR = factor(CHR, levels=new))

ggplot(plot_af, aes(x = bin / 1e6, y = value)) +
  geom_line(linewidth = 0.4, colour = "grey40") +
  geom_point(size = 0.6, colour = "steelblue") +
  facet_grid(rows = vars(CHR)) +
  labs(
    x = "Position (Mb)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text.y = element_text(size = 10),
    axis.text.y = element_text(size=8),
    panel.spacing.y = unit(0.5, "lines"),
    strip.background = element_rect(color="black", size=1, linetype="solid"),
    strip.text = element_text(size=15),
    axis.title.x = element_text(size=13)
  )

# Remove extranuclear DNA
df_exnuc_rem = filter(plot_dt, CHR != "HG994977.1" & CHR != "HG994976.1")

df_exnuc_rem$CHR = factor(df_exnuc_rem$CHR, levels = unique(df_exnuc_rem$CHR))

ggplot(df_exnuc_rem, aes(x = bin / 1e6, y = value)) +
  geom_line(size=0.4) +
  facet_wrap(~ CHR,
             scales = "free_x",
             nrow = 3) +
  labs(
    x = "\nGenomic position (Mb)",
    y = "Mean \u0394 SNP Index\n"
  ) +
  geom_hline(yintercept=0.9, 
             linetype="dashed", 
             color = "red",
             size=0.5) +
  scale_y_continuous(breaks = sort(c(seq(0, 1, length.out=5), 0.9))) +
  theme_bw() + 
  theme(strip.background = element_rect(colour="black", fill="#dddddd",size=1.2),
        strip.text = element_text(size = 15),
        axis.title.y = element_text(margin = margin(t = 20,b = 0), size = 20),
        axis.title.x = element_text(margin = margin(t = 10,b = 0), size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
