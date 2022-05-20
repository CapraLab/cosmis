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

# make a density plot
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme_classic() + theme(
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
  plot.margin = plot_margin
)


# make a violin plot
fig_s6a <- combined %>% 
  ggplot(
    mapping = aes(
      x = factor(str_source, levels = c("PDB", "SWISS-MODEL", "AF2")), 
      y = cosmis, 
      fill = factor(str_source, levels = c("PDB", "SWISS-MODEL", "AF2"))
    )
  ) +
  geom_violin(
    trim = FALSE,
    size = 1,
    alpha = 0.5
  ) +
  geom_boxplot(
    width = 0.2,
    fill = "white",
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = pnw_palette(name = "Shuksan2", n = 7)[c(2, 4, 6)]
  ) + 
  labs(
    y = "COSMIS score"
  ) +
  scale_x_discrete(
    labels = c("PDB", "SWISS-MODEL", "AF2")
  ) +
  geom_abline(
    slope = 0,
    intercept = 0, 
    size = 0.5,
    colour = "gray",
    linetype = "dashed"
  ) + 
  scale_y_continuous(
    breaks = seq(-8, 8, 2),
    limits = c(-8, 8, 2),
    labels = sprintf("%.1f", seq(-8, 8, 2)),
    expand = c(0, 0)
  ) +
  plot_theme + theme(
    axis.title.x = element_blank()
  )

ggsave(
  filename = paste(data_dir, "fig_s6", "fig_s6a.svg", sep = "/"),
  plot = fig_s6a,
  width = 8,
  height = 6,
  units = "in",
  device = "svg",
)

# make a violin plot
fig_s6b <- combined %>% 
  ggplot(
    mapping = aes(
      x = factor(str_source, levels = c("PDB", "SWISS-MODEL", "AF2")), 
      y = cossyn, 
      fill = factor(str_source, levels = c("PDB", "SWISS-MODEL", "AF2"))
    )
  ) +
  geom_violin(
    trim = FALSE,
    size = 1,
    alpha = 0.5
  ) +
  geom_boxplot(
    width = 0.2,
    fill = "white",
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = pnw_palette(name = "Shuksan2", n = 7)[c(2, 4, 6)]
  ) + 
  labs(
    y = "Contact set synonymous tolerance score"
  ) +
  scale_x_discrete(
    labels = c("PDB", "SWISS-MODEL", "AF2")
  ) +
  geom_abline(
    slope = 0,
    intercept = 0, 
    size = 0.5,
    colour = "gray",
    linetype = "dashed"
  ) + 
  scale_y_continuous(
    breaks = seq(-8, 8, 2),
    limits = c(-8, 8, 2),
    labels = sprintf("%.1f", seq(-8, 8, 2)),
    expand = c(0, 0)
  ) +
  plot_theme + theme(
    axis.title.x = element_blank()
  )

ggsave(
  filename = paste(data_dir, "fig_s6", "fig_s6b.svg", sep = "/"),
  plot = fig_s6b,
  width = 8,
  height = 6,
  units = "in",
  device = "svg",
)

