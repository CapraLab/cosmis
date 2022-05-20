# load required libraries
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

# compute residue-level coverages
pdb_cov <- cosmis_df_pdb %>% group_by(uniprot_id) %>% summarise(
  coverage = n() / max(uniprot_length)
)
sm_cov <- cosmis_df_swiss_model %>% group_by(uniprot_id) %>% summarise(
  coverage = n() / max(uniprot_length)
)
af2_cov <- cosmis_df_alphafold %>% group_by(uniprot_id) %>% summarise(
  coverage = n() / max(uniprot_length)
)

cov_tib <- bind_rows(
  tibble(
    cov = pdb_cov$coverage,
    dataset = "PDB"
  ),
  tibble(
    cov = sm_cov$coverage,
    dataset = "SWISS-MODEL"
  ),
  tibble(
    cov = af2_cov$coverage,
    dataset = "AF2"
  )
)

# make a density plot
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme_classic() + theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text = element_text(
    size = 20, color = "black"
  ),
  axis.text.x = element_text(
    margin = margin(t = 10)
  ),
  axis.title.x = element_text(
    color = "black", size = 24, margin = margin(t = 10)
  ),
  axis.title.y = element_text(
    color = "black", size = 24, margin = margin(r = 10)
  ),
  axis.ticks.length = unit(.25, "cm"),
  axis.ticks = element_line(),
  legend.position = "none",
  plot.margin = plot_margin
)

# fraction of transcript residues covered by structure
fig_s1a <- cov_tib %>% 
  ggplot(
    mapping = aes(
      x = factor(dataset, levels = c("PDB", "SWISS-MODEL", "AF2")), 
      y = cov, 
      fill = factor(dataset, levels = c("PDB", "SWISS-MODEL", "AF2"))
    )
  ) +
  geom_violin(
    trim = FALSE,
    size = 1
  ) +
  geom_boxplot(
    width = 0.1,
    fill = "white",
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = pnw_palette(name = "Shuksan2", n = 7)[c(2, 4, 6)]
  ) + 
  labs(
    x = "",
    y = str_wrap("Fraction of residues covered by structure", width = 30)
  ) +
  scale_x_discrete(
    labels = c("PDB", "SWISS-MODEL", "AF2")
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1), 
    labels = sprintf("%.1f", seq(0, 1, 0.2)),
    expand = c(0, 0)
  ) +
  plot_theme

# save the plot
ggsave(
  filename = paste(data_dir, "fig_s1", "fig_s1a.svg", sep = "/"),
  plot = fig_s1a,
  width = 8,
  height = 6,
  units = "in",
  device = "svg",
)
