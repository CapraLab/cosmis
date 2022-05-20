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

#===============================================================================
# cumulative distribution of sequence separation
#===============================================================================

# get column of sequence distance
get_sep_col <- function(data, low = 5, high = Inf) {
  x <- data %>% pull(seq_separations) %>% 
    str_split(";") %>% 
    unlist() %>% 
    parse_integer() %>% 
    enframe(name = "index", value = "sep") %>% 
    # each contact was counted twice, once for A->B and once for B->A
    # retain only one instance
    filter(sep > low & sep <= high)
  return(x)
}

pdb_seq_sep <- get_sep_col(cosmis_df_pdb, low = 0)
sm_seq_sep <- get_sep_col(cosmis_df_swiss_model, low = 0)
af2_seq_sep <- get_sep_col(cosmis_df_alphafold, low = 0)

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
  aspect.ratio = 1.0,
  plot.margin = plot_margin
)

# PDB
fig_s3a_pdb <- ggplot() +
  stat_ecdf(
    data = pdb_seq_sep,
    mapping = aes(x = sep),
    geom = "smooth",
    color = pnw_palette(name = "Shuksan2", n = 7)[2],
    size = 2
  ) +
  geom_vline(
    xintercept = 15, 
    color = "blue", 
    linetype = "dotted", 
    size = 1
  ) +
  labs(
    x = "Separation in 1D sequence (# residues)", 
    y = "Cumulative fraction of contacts"
  ) +
  scale_x_continuous(
    limits = c(0, 500),
    breaks = seq(0, 500, 100),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1.2),
    breaks = seq(0, 1.2, 0.2)
  ) + 
  plot_theme

ggsave(
  file = str_c(data_dir, "fig_s3", "fig_s3a_pdb.svg", sep = "/"),
  plot = fig_s3a_pdb,
  device = "svg",
  unit = "in",
  width = 6,
  height = 6
)

# SWISS-MODEL
fig_s3a_swiss_model <- ggplot() +
  stat_ecdf(
    data = sm_seq_sep,
    mapping = aes(x = sep),
    geom = "smooth",
    color = pnw_palette(name = "Shuksan2", n = 7)[4],
    size = 2
  ) +
  geom_vline(
    xintercept = 15, 
    color = "blue", 
    linetype = "dotted", 
    size = 1
  ) +
  labs(
    x = "Separation in 1D sequence (# residues)", 
    y = "Cumulative fraction of contacts"
  ) +
  scale_x_continuous(
    limits = c(0, 500),
    breaks = seq(0, 500, 100),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1.2),
    breaks = seq(0, 1.2, 0.2)
  ) + 
  plot_theme

ggsave(
  file = str_c(data_dir, "fig_s3", "fig_s3a_swiss_model.svg", sep = "/"),
  plot = fig_s3a_swiss_model,
  device = "svg",
  unit = "in",
  width = 6,
  height = 6
)

# AlphaFold
fig_s3a_alphafold <- ggplot() +
  stat_ecdf(
    data = af2_seq_sep,
    mapping = aes(x = sep),
    geom = "smooth",
    color = pnw_palette(name = "Shuksan2", n = 7)[6],
    size = 2
  ) +
  geom_vline(
    xintercept = 15, 
    color = "blue", 
    linetype = "dotted", 
    size = 1
  ) +
  labs(
    x = "Separation in 1D sequence (# residues)", 
    y = "Cumulative fraction of contacts"
  ) +
  scale_x_continuous(
    limits = c(0, 500),
    breaks = seq(0, 500, 100),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1.2),
    breaks = seq(0, 1.2, 0.2)
  ) + 
  plot_theme

ggsave(
  file = str_c(data_dir, "fig_s3", "fig_s3a_alphafold.svg", sep = "/"),
  plot = fig_s3a_alphafold,
  device = "svg",
  unit = "in",
  width = 6,
  height = 6
)