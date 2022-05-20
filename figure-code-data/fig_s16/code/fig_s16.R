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

# combine the data set
comm_cols <- Reduce(
  intersect,
  list(
    colnames(cosmis_df_pdb), 
    colnames(cosmis_df_swiss_model),
    colnames(cosmis_df_alphafold)
  )
)
combined <- rbind(
  cosmis_df_pdb[comm_cols],
  cosmis_df_swiss_model[comm_cols],
  cosmis_df_alphafold[comm_cols]
)

plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
hist_plot_theme <- theme_classic() + theme(
  panel.border = element_blank(),
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
  axis.ticks = element_line(),
  legend.position = "none",
  plot.margin = plot_margin
)

fig_s16a <-  ggplot(data = combined) + 
  geom_histogram(
    mapping = aes(cs_syn_obs / cs_syn_poss),
    binwidth = 0.025,
    boundary = 0,
    colour = "black",
    size = 0.3,
    fill = pnw_palette(name = "Shuksan2", n = 7)[2]
  ) +
  geom_vline(
    xintercept = 0, color = "grey", linetype = "dashed"
  ) +
  scale_x_continuous(
    name = str_wrap(
      "Proportion of possible synonymous variants", width = 30
    ),
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = str_wrap("Number of contact sets (million)", width = 30),
    limits = c(0, 2000000),
    breaks = seq(0, 2000000, 500000),
    labels = c("0", "0.5", "1.0", "1.5", "2.0"),
    expand = c(0, 0)
  ) +
  hist_plot_theme

ggsave(
  filename = paste(data_dir, "fig_s16", "fig_s16a.svg", sep = "/"),
  plot = fig_s16a,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)

fig_s16b <-  ggplot(data = combined) + 
  geom_histogram(
    mapping = aes(cs_mis_obs / cs_mis_poss),
    binwidth = 0.025,
    boundary = 0,
    colour = "black",
    size = 0.3,
    fill = pnw_palette(name = "Shuksan2", n = 7)[6]
  ) +
  geom_vline(
    xintercept = 0, color = "grey", linetype = "dashed"
  ) +
  scale_x_continuous(
    name = str_wrap(
      "Proportion of possible missense variants", width = 30
    ),
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = str_wrap("Number of contact sets (million)", width = 30),
    limits = c(0, 2000000),
    breaks = seq(0, 2000000, 500000),
    labels = c("0", "0.5", "1.0", "1.5", "2.0"),
    expand = c(0, 0)
  ) +
  hist_plot_theme

ggsave(
  filename = paste(data_dir, "fig_s16", "fig_s16b.svg", sep = "/"),
  plot = fig_s16b,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)
