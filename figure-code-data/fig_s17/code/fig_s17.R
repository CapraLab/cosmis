# load required packages
rm(list = ls())
library(biomaRt)
library(ggplot2)
library(tidyverse)
library(PNWColors)

# path to source data folder
data_dir <- "/path/to/source_data"

prob_vs_count <- read.table(
  file = paste(
    data_dir, "fig_s17/data", 
    "gnomad_all_enst_obs_prob_vs_counts.tsv", 
    sep = "/"
  ),
  header = TRUE
)

samocha_rates <- read.csv(
  file = paste(
    data_dir, "fig_s17/data",
    "samocha_mutation_rates.csv", sep = "/"
  ),
  header = TRUE,
  row.names = "gene"
)

# map ensembl transcript identifiers to refseq mRNA identifers
ensembl <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
)
map <- getBM(
  attributes = c("ensembl_transcript_id", "refseq_mrna", "hgnc_symbol"),
  filters = c("ensembl_transcript_id"),
  values = prob_vs_count$enst_id,
  mart = ensembl
)
map[map == ""] <- NA
map <- na.omit(map)

# ensembl transcript to hgnc_symbol mapping, duplicates removed
enst_to_hgnc <- map[, c("ensembl_transcript_id", "hgnc_symbol")]
enst_to_hgnc <- enst_to_hgnc[!duplicated(enst_to_hgnc), ]
common_genes <- enst_to_hgnc[
  enst_to_hgnc$hgnc_symbol %in% rownames(samocha_rates), 
]

# get mutation rates from samocha
cnl_enst_to_hgnc <- read_tsv(
  file = paste(
    data_dir, "fig_s17/data", 
    "cnl_enst_to_hgnc.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)
samocha_syn_rates <- samocha_rates[cnl_enst_to_hgnc$hgnc_symbol, "syn"]
samocha_mis_rates <- samocha_rates[cnl_enst_to_hgnc$hgnc_symbol, "mis"]

rownames(prob_vs_count) <- prob_vs_count$enst_id
my_syn_rates <- prob_vs_count[
  cnl_enst_to_hgnc$ensembl_transcript_id,
  "syn_prob"
]
my_mis_rates <- prob_vs_count[
  cnl_enst_to_hgnc$ensembl_transcript_id,
  "mis_prob"
]

mutation_rates <- tibble(
  my_syn_rates = log10(my_syn_rates),
  my_mis_rates = log10(my_mis_rates),
  samocha_syn_rates,
  samocha_mis_rates
) %>% drop_na()

# make a scatter plot
plot_margin <- margin(
  t = 0.5, r = 1.0, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text = element_text(
    size = 16, color = "black"
  ),
  axis.title.x = element_text(
    color = "black", size = 20, margin = margin(t = 10)
  ),
  axis.title.y = element_text(
    color = "black", size = 20, margin = margin(r = 10)
  ),
  axis.ticks.length = unit(.25, "cm"),
  axis.ticks = element_line(),
  legend.position = "none",
  aspect.ratio = 1.0,
  plot.margin = plot_margin
)

fig_s17a <- ggplot(
  data = mutation_rates,
  mapping = aes(
    x = my_syn_rates,
    y = samocha_syn_rates
  )
) +
  geom_point(
    pch = 19,
    size = 1.5,
    color = pnw_palette(name = "Shuksan2", n = 7)[2],
    alpha = 0.5
  ) +
  scale_x_continuous(
    name = "",
    limits = c(-8, -2),
    breaks = seq(-8, -2, 2),
    labels = c("-8.0", "-6.0", "-4.0", "-2.0"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "",
    limits = c(-8, -2),
    breaks = seq(-8, -2, 2),
    labels = c("-8.0", "-6.0", "-4.0", "-2.0"),
    expand = c(0, 0)
  ) +
  plot_theme

# save the scatter plot to disk
ggsave(
  filename = paste(data_dir, "fig_s17", "fig_s17a.svg", sep = "/"),
  plot = fig_s17a,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)

fig_s17b <- ggplot(
  data = mutation_rates,
  mapping = aes(
    x = my_mis_rates,
    y = samocha_mis_rates
  )
) +
  geom_point(
    pch = 19,
    size = 1.5,
    color = pnw_palette(name = "Shuksan2", n = 7)[6],
    alpha = 0.5
  ) +
  scale_x_continuous(
    name = "",
    limits = c(-8, -2),
    breaks = seq(-8, -2, 2),
    labels = c("-8.0", "-6.0", "-4.0", "-2.0"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "",
    limits = c(-8, -2),
    breaks = seq(-8, -2, 2),
    labels = c("-8.0", "-6.0", "-4.0", "-2.0"),
    expand = c(0, 0)
  ) +
  plot_theme

# save the scatter plot to disk
ggsave(
  filename = paste(data_dir, "fig_s17", "fig_s17b.svg", sep = "/"),
  plot = fig_s17b,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)