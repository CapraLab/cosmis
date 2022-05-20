# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())


# path to source data folder
data_dir <- "/path/to/source_data/"

# load data sets
cosmis_pdb <- read_tsv(
  file = paste(data_dir, "cosmis_scores_pdb.tsv.gz", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    cosmis = col_double(),
    cosmis_pvalue = col_double(),
    cossyn = col_double(),
    str_source = col_character()
  )
)
cosmis_swiss_model <- read_tsv(
  file = paste(data_dir, "cosmis_scores_swiss_model.tsv.gz", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    cosmis = col_double(),
    cosmis_pvalue = col_double(),
    cossyn = col_double(),
    str_source = col_character()
  )
)
cosmis_alphafold <- read_tsv(
  file = paste(data_dir, "cosmis_scores_alphafold.tsv.gz", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    cosmis = col_double(),
    cosmis_pvalue = col_double(),
    cossyn = col_double(),
    str_source = col_character()
  )
)
combined <- rbind(
  cosmis_pdb,
  cosmis_swiss_model,
  cosmis_alphafold
)


# setting plot parameters
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
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

fig_3b <- combined %>% ggplot() +
  geom_density(
    mapping = aes(
      x = cosmis
    ),
    fill = pnw_palette(name = "Shuksan2", n = 7)[7],
    color = pnw_palette(name = "Shuksan2", n = 7)[7],
    alpha = 0.5,
    size = 1.2
  ) +
  scale_x_continuous(
    name = str_wrap("COSMIS score", width = 30),
    limits = c(-8, 8),
    breaks = seq(-8, 8, 4),
    labels = sprintf("%.1f", seq(-8, 8, 4)),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Density",
    limits = c(0, 0.4),
    breaks = seq(0, 0.4, 0.1),
    labels = sprintf("%.2f", seq(0, 0.4, 0.1)),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_3b", "fig_3b.svg", sep = "/"),
  plot = fig_3b,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)