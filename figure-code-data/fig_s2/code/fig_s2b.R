# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())

# path to source data folder
data_dir <- "/path/to/source_data"

# load data sets
cosmis_df_pdb <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_pdb.tsv.gz", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    pdb_pos = col_integer(),
    pdb_aa = col_character(),
    pdb_id = col_character(),
    chain_id = col_character(),
    seq_separations = col_character(),
    num_contacts = col_integer(),
    cs_syn_poss = col_integer(),
    cs_mis_poss = col_integer(),
    cs_gc_content = col_double(),
    cs_syn_prob = col_double(),
    cs_syn_obs = col_integer(),
    cs_mis_prob = col_double(),
    cs_mis_obs = col_integer(),
    mis_pmt_mean = col_double(),
    mis_pmt_sd = col_double(),
    mis_p_value = col_double(),
    syn_pmt_mean = col_double(),
    syn_pmt_sd = col_double(),
    syn_p_value = col_double(),
    enst_syn_obs = col_double(),
    enst_mis_obs = col_double(),
    enst_syn_exp = col_double(),
    enst_mis_exp = col_double(),
    uniprot_length = col_integer()
  )
)

# load swiss model dataset
cosmis_df_swiss_model <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_swiss_model.tsv.gz", sep = "/"),
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    seq_separations = col_character(),
    num_contacts = col_integer(),
    cs_syn_poss = col_integer(),
    cs_mis_poss = col_integer(),
    cs_gc_content = col_double(),
    cs_syn_prob = col_double(),
    cs_syn_obs = col_integer(),
    cs_mis_prob = col_double(),
    cs_mis_obs = col_integer(),
    mis_pmt_mean = col_double(),
    mis_pmt_sd = col_double(),
    mis_p_value = col_double(),
    syn_pmt_mean = col_double(),
    syn_pmt_sd = col_double(),
    syn_p_value = col_double(),
    enst_syn_obs = col_double(),
    enst_mis_obs = col_double(),
    enst_syn_exp = col_double(),
    enst_mis_exp = col_double(),
    uniprot_length = col_integer()
  )
)

# load alphafold dataset
cosmis_df_alphafold <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_alphafold.tsv.gz", sep = "/"),
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    seq_separations = col_character(),
    num_contacts = col_integer(),
    cs_syn_poss = col_integer(),
    cs_mis_poss = col_integer(),
    cs_gc_content = col_double(),
    cs_syn_prob = col_double(),
    cs_syn_obs = col_integer(),
    cs_mis_prob = col_double(),
    cs_mis_obs = col_integer(),
    mis_pmt_mean = col_double(),
    mis_pmt_sd = col_double(),
    mis_p_value = col_double(),
    syn_pmt_mean = col_double(),
    syn_pmt_sd = col_double(),
    syn_p_value = col_double(),
    enst_syn_obs = col_double(),
    enst_mis_obs = col_double(),
    enst_syn_exp = col_double(),
    enst_mis_exp = col_double(),
    uniprot_length = col_integer()
  )
)

# variant count statistics by contact set size
pdb_by_cs_size <- cosmis_df_pdb %>% group_by(num_contacts) %>% 
  summarise(
    mean_syn = mean(cs_syn_obs),
    sd_syn = sd(cs_syn_obs),
    mean_mis = mean(cs_mis_obs),
    sd_mis = sd(cs_mis_obs),
    count = n()
  )

sm_by_cs_size <- cosmis_df_swiss_model %>% group_by(num_contacts) %>% 
  summarise(
    mean_syn = mean(cs_syn_obs),
    sd_syn = sd(cs_syn_obs),
    mean_mis = mean(cs_mis_obs),
    sd_mis = sd(cs_mis_obs),
    count = n()
  )

af2_by_cs_size <- cosmis_df_alphafold %>% group_by(num_contacts) %>% 
  summarise(
    mean_syn = mean(cs_syn_obs),
    sd_syn = sd(cs_syn_obs),
    mean_mis = mean(cs_mis_obs),
    sd_mis = sd(cs_mis_obs),
    count = n()
  )

# plot theme
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

fig_s2b_pdb <- pdb_by_cs_size %>% 
  filter(
    num_contacts != 2, num_contacts < 22
  ) %>% ggplot() +
  geom_pointrange(
    mapping = aes(
      x = num_contacts, y = mean_syn,
      ymin = mean_syn - sd_syn, ymax = mean_syn + sd_syn),
      color = pnw_palette(name = "Shuksan2", n = 5)[1]
  ) +
  geom_pointrange(
    mapping = aes(
      x = num_contacts, y = mean_mis,
      ymin = mean_mis - sd_mis, ymax = mean_mis + sd_mis),
      color = pnw_palette(name = "Shuksan2", n = 5)[5]
  ) +
  scale_x_continuous(
    name = "Contact set size",
    limits = c(0, 25),
    breaks = seq(0, 25, 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Variant count",
    limits = c(0, 16),
    breaks = seq(0, 16, 4),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s2", "fig_s2b_pdb.svg", sep = "/"),
  plot = fig_s2b_pdb,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)

fig_s2b_swiss_model <- sm_by_cs_size %>% 
  filter(
    num_contacts != 2, num_contacts < 22
  ) %>% ggplot() +
  geom_pointrange(
    mapping = aes(
      x = num_contacts, y = mean_syn,
      ymin = mean_syn - sd_syn, ymax = mean_syn + sd_syn),
    color = pnw_palette(name = "Shuksan2", n = 5)[1]
  ) +
  geom_pointrange(
    mapping = aes(
      x = num_contacts, y = mean_mis,
      ymin = mean_mis - sd_mis, ymax = mean_mis + sd_mis),
    color = pnw_palette(name = "Shuksan2", n = 5)[5]
  ) +
  scale_x_continuous(
    name = "Contact set size",
    limits = c(0, 25),
    breaks = seq(0, 25, 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Variant count",
    limits = c(0, 16),
    breaks = seq(0, 16, 4),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s2", "fig_s2b_swiss_model.svg", sep = "/"),
  plot = fig_s2b_swiss_model,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)

fig_s2b_alphafold <- af2_by_cs_size %>% 
  filter(
    num_contacts > 2, num_contacts < 22
  ) %>% ggplot() +
  geom_pointrange(
    mapping = aes(
      x = num_contacts, y = mean_syn,
      ymin = mean_syn - sd_syn, ymax = mean_syn + sd_syn),
    color = pnw_palette(name = "Shuksan2", n = 5)[1]
  ) +
  geom_pointrange(
    mapping = aes(
      x = num_contacts, y = mean_mis,
      ymin = mean_mis - sd_mis, ymax = mean_mis + sd_mis),
    color = pnw_palette(name = "Shuksan2", n = 5)[5]
  ) +
  scale_x_continuous(
    name = "Contact set size",
    limits = c(0, 25),
    breaks = seq(0, 25, 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Variant count",
    limits = c(0, 16),
    breaks = seq(0, 16, 4),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s2", "fig_s2b_alphafold.svg", sep = "/"),
  plot = fig_s2b_alphafold,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)
