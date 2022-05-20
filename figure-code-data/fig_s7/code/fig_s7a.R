# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())

# path to source data folder
data_dir <- "/path/to/source_data"

pdb <- read_tsv(
  file = str_c(
    data_dir, "fig_s7/data", "pdb_hc_sites_cluster.tsv", sep = "/"
  ),
  col_names = TRUE,
  na = "nan"
)
swiss_model <- read_tsv(
  file = str_c(
    data_dir, "fig_s7/data", "sm_hc_sites_cluster.tsv", sep = "/"
  ),
  col_names = TRUE,
  na = "nan"
)
alphafold <- read_tsv(
  file = str_c(
    data_dir, "fig_s7/data", "af2_hc_sites_cluster.tsv", sep = "/"
  ), 
  col_names = TRUE,
  na = "nan"
)

# plot theme
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme_classic() + theme(
  panel.border = element_rect(
    colour = "black", 
    size = 1, 
    fill = "transparent"
  ),
  axis.text = element_text(
    size = 16, color = "black"
  ),
  axis.title.x = element_text(
    color = "black", size = 20, margin = margin(t = 10, r = 10)
  ),
  axis.title.y = element_text(
    color = "black", size = 20, margin = margin(r = 10)
  ),
  axis.ticks.length = unit(.25, "cm"),
  axis.ticks = element_line(size = 1.0),
  legend.position = "none",
  aspect.ratio = 1.0,
  plot.margin = plot_margin
)

# PDB
fig_s7a_pdb <- pdb %>% filter(n_sites >= 2) %>% 
  mutate(
    is_sig = ifelse(p_value < 0.01, "YES", "NO")
  ) %>% 
  ggplot() + 
  geom_point(
    aes(x = mean_d, y = permuted_d, color = is_sig), size = 1
  ) + 
  geom_abline(
    intercept = 0, slope = 1
  ) + 
  scale_x_continuous(
    name = str_wrap("Avg. pairwise distance (Å)"),
    breaks = seq(0, 125, 25),
    limits = c(0, 125),
    expand = c(0, 0)
  ) + 
  scale_y_continuous(
    name = str_wrap("Avg. pairwise distance (Å) (1000 permutations)", width = 30),
    breaks = seq(0, 125, 25),
    limits = c(0, 125),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = pnw_palette(name = "Shuksan2", n = 7)[c(2, 6)]
  ) + 
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s7", "fig_s7a_pdb.svg", sep = "/"),
  plot = fig_s7a_pdb,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)

# SWISS-MODEL
fig_s7a_swiss_model <- swiss_model %>% filter(n_sites >= 2) %>% 
  mutate(
    is_sig = ifelse(p_value < 0.01, "YES", "NO")
  ) %>% 
  ggplot() + 
  geom_point(
    aes(x = mean_d, y = permuted_d, color = is_sig), size = 1
  ) + 
  geom_abline(
    intercept = 0, slope = 1
  ) + 
  scale_x_continuous(
    name = str_wrap("Avg. pairwise distance (Å)"),
    breaks = seq(0, 125, 25),
    limits = c(0, 125),
    expand = c(0, 0)
  ) + 
  scale_y_continuous(
    name = str_wrap("Avg. pairwise distance (Å) (1000 permutations)", width = 30),
    breaks = seq(0, 125, 25),
    limits = c(0, 125),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = pnw_palette(name = "Shuksan2", n = 7)[c(2, 6)]
  ) + 
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s7", "fig_s7a_swiss_model.svg", sep = "/"),
  plot = fig_s7a_swiss_model,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)

# AlphaFold
fig_s7a_alphafold <- alphafold %>% filter(n_sites >= 2) %>% 
  mutate(
    is_sig = ifelse(p_value < 0.01, "YES", "NO")
  ) %>% 
  ggplot() + 
  geom_point(
    aes(x = mean_d, y = permuted_d, color = is_sig), size = 1
  ) + 
  geom_abline(
    intercept = 0, slope = 1
  ) + 
  scale_x_continuous(
    name = str_wrap("Avg. pairwise distance (Å)"),
    breaks = seq(0, 125, 25),
    limits = c(0, 125),
    expand = c(0, 0)
  ) + 
  scale_y_continuous(
    name = str_wrap("Avg. pairwise distance (Å) (1000 permutations)", width = 30),
    breaks = seq(0, 125, 25),
    limits = c(0, 125),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = pnw_palette(name = "Shuksan2", n = 7)[c(2, 6)]
  ) + 
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s7", "fig_s7a_alphafold.svg", sep = "/"),
  plot = fig_s7a_alphafold,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)