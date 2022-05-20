# load required packages
rm(list = ls())
library(tidyverse)
library(PNWColors)

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

pdb_cs_size <- data.frame(
  size = cosmis_df_pdb$num_contacts, dataset = "PDB"
)
sm_cs_size <- data.frame(
  size = cosmis_df_swiss_model$num_contacts, dataset = "SWISS-MODEL"
)
af2_cs_size <- data.frame(
  size = cosmis_df_alphafold$num_contacts, dataset = "AF2"
)
cs_size = rbind(pdb_cs_size, sm_cs_size, af2_cs_size)

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


fig_s2a_pdb <- ggplot(
    data = pdb_cs_size,
    mapping = aes(x = size, fill = dataset)
  ) +
  geom_bar(
    color = "black",
    alpha = 0.5,
    position = "dodge"
  ) +
  scale_fill_manual(
      values = pnw_palette(name = "Shuksan2", n = 7)[c(2)]
  ) +
  scale_x_continuous(
    name = "Contact set size",
    limits = c(0, 25),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = expression(paste("Frequency / 1000")),
    limits = c(0, 200000),
    labels = seq(0, 200, 50),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s2", "fig_s2a_pdb.svg", sep = "/"),
  plot = fig_s2a_pdb,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)


fig_s2a_swiss_model <- ggplot(
    data = sm_cs_size,
    mapping = aes(x = size, fill = dataset)
  ) +
  geom_bar(
    color = "black",
    alpha = 0.5,
    position = "dodge"
  ) +
  scale_fill_manual(
    values = pnw_palette(name = "Shuksan2", n = 7)[c(4)]
  ) +
  scale_x_continuous(
    name = "Contact set size",
    limits = c(0, 25),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = expression(paste("Frequency / 1000")),
    limits = c(0, 250000),
    labels = seq(0, 250, 50),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s2", "fig_s2a_swiss_model.svg", sep = "/"),
  plot = fig_s2a_swiss_model,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)


fig_s2a_alphafold <- ggplot(
  data = af2_cs_size,
  mapping = aes(x = size, fill = dataset)
) +
  geom_bar(
    color = "black",
    alpha = 0.5,
    position = "dodge"
  ) +
  scale_fill_manual(
    values = pnw_palette(name = "Shuksan2", n = 7)[c(6)]
  ) +
  scale_x_continuous(
    name = "Contact set size",
    limits = c(0, 25),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = expression(paste("Frequency / 1000")),
    limits = c(0, 800000),
    breaks = seq(0, 800000, 100000),
    labels = seq(0, 800, 100),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s2", "fig_s2a_alphafold.svg", sep = "/"),
  plot = fig_s2a_alphafold,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)
