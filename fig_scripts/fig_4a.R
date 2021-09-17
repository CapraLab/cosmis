# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())

# data and figure paths
data_dir <- "/path/to/datasets/"
figure_dir <- "/path/to/where/figures/are/stored/"

# load data sets
cosmis_df_pdb <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_pdb.tsv", sep = "/"), 
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
    syn_var_sites = col_integer(),
    total_syn_sites = col_double(),
    mis_var_sites = col_integer(),
    total_mis_sites = col_double(),
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
  file = paste(data_dir, "cosmis_dataset_swiss_model.tsv", sep = "/"),
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    seq_separations = col_character(),
    num_contacts = col_integer(),
    syn_var_sites = col_integer(),
    total_syn_sites = col_double(),
    mis_var_sites = col_integer(),
    total_mis_sites = col_double(),
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
  file = paste(data_dir, "cosmis_dataset_alphafold.tsv", sep = "/"),
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    seq_separations = col_character(),
    num_contacts = col_integer(),
    syn_var_sites = col_integer(),
    total_syn_sites = col_double(),
    mis_var_sites = col_integer(),
    total_mis_sites = col_double(),
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


# load pLI scores
pli_pdb <- read_csv(
  file = paste(data_dir, "pdb_dataset_pli_mis_z.csv", sep = "/"),
  col_names = TRUE
)
pli_sm <- read_csv(
  file = paste(data_dir, "sm_dataset_pli_mis_z.csv", sep = "/"),
  col_names = TRUE
)
pli_af2 <- read_csv(
  file = paste(data_dir, "af2_dataset_pli_mis_z.csv", sep = "/"),
  col_names = TRUE
)
# combine pLI scores
pli_combined <- bind_rows(
  pli_pdb, pli_sm, pli_af2
)

# combine pLI and COSMIS datasets
cosmis_pli_combined <- bind_cols(
  combined, pli_combined
)

cosmis_vs_pli <- cosmis_pli_combined %>% drop_na %>% 
  group_by(uniprot_id) %>% summarise(
    n_res = n(),
    mean_cosmis = mean(z_score),
    mean_pli = mean(pli)
  ) %>%  mutate(class = ifelse(
    mean_pli >= 0.9,
    "Intolerant",
    ifelse(mean_pli <= 0.1, "Tolerant", "Unsure")
  )
  )

# make violin plot
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
  plot.margin = plot_margin
)

cosmis_vs_pli_plot <- cosmis_vs_pli %>% 
  ggplot(
    aes(x = reorder(class, mean_cosmis), y = mean_cosmis, fill = class)
  ) + 
  geom_violin(
    size = 0.8
  ) +
  scale_fill_manual(
    values = pnw_palette(name = "Shuksan2", n = 5)[c(4, 2, 3)]
  ) +
  geom_boxplot(
    size = 0.8,
    width = 0.12,
    fill = "white",
    outlier.shape = NA
  ) +
  scale_x_discrete(
    name = ""
  ) +
  scale_y_continuous(
    name = "Mean COSMIS score",
    limits = c(-4, 4),
    labels = sprintf("%.1f", seq(-4, 4, 2)),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(figure_dir, "cosmis_vs_pli_plot.svg", sep = "/"),
  plot = cosmis_vs_pli_plot,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)