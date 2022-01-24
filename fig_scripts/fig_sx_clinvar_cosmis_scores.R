library(tidyverse)
library(PNWColors)
rm(list = ls())

# data and figure paths
# data_dir <- "/path/to/datasets"
# figure_dir <- "/path/to/where/figures/are/stored"

data_dir <- "/Users/lib14/OneDrive/Research/projects/cosmis/results"
figure_dir <- "/Users/lib14/OneDrive/Research/projects/cosmis/figures"

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
combined <- combined %>% mutate(
  cosmis = (cs_mis_obs - mis_pmt_mean) / mis_pmt_sd,
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)

# load variant IDs
clinvar_patho_ids <-  read_tsv(
  file = paste(
    data_dir, 
    "clinvar_unambiguous_pathogenic_ids_20210807.tsv", 
    sep = "/"
  ), 
  col_names = TRUE
)
clinvar_patho_ids <- clinvar_patho_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)
clinvar_benign_ids <- read_tsv(
  file = paste(
    data_dir, 
    "clinvar_unambiguous_benign_ids_20210807.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)
clinvar_benign_ids <- clinvar_benign_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)
clinvar_likely_patho_ids <-  read_tsv(
  file = paste(
    data_dir, 
    "clinvar_unambiguous_likely_pathogenic_ids_20210807.tsv", 
    sep = "/"
  ), 
  col_names = TRUE
)
clinvar_likely_patho_ids <- clinvar_likely_patho_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)
clinvar_likely_benign_ids <- read_tsv(
  file = paste(
    data_dir, 
    "clinvar_unambiguous_likely_benign_ids_20210807.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)
clinvar_likely_benign_ids <- clinvar_likely_benign_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)

# cosmis scores for ClinVar pathogenic variants
patho_cosmis <- combined %>% filter(
  full_id %in% clinvar_patho_ids$full_id
) %>% mutate(
  class = "pathogenic"
)
benign_cosmis <- combined %>% filter(
  full_id %in% clinvar_benign_ids$full_id
) %>% mutate(
  class = "benign"
)
likely_patho_cosmis <- combined %>% filter(
  full_id %in% clinvar_likely_patho_ids$full_id
) %>% mutate(
  class = "likely_pathogenic"
)
likely_benign_cosmis <- combined %>% filter(
  full_id %in% clinvar_likely_benign_ids$full_id
) %>% mutate(
  class = "likely benign"
)

cosmis_clinvar <- rbind(
  patho_cosmis,
  benign_cosmis,
  likely_patho_cosmis,
  likely_benign_cosmis
)
cosmis_clinvar <- cosmis_clinvar %>% mutate(
  class = as_factor(class)
)

# wilcox test
# patho <- subset(x = cosmis_clinvar, subset = label == 1)
# benign <- subset(x = cosmis_clinvar, subset = label == 0)
w_test <- wilcox.test(likely_patho_cosmis$cosmis, patho_cosmis$cosmis)

# set plot theme
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme_classic() + theme(
  axis.text.x = element_text(
    size = 16, color = "black", angle = 90, hjust = 1, vjust = 0.5
  ),
  axis.text.y = element_text(
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
cosmis_violin_clinvar <- cosmis_clinvar %>% 
  filter(class %in% c("benign", "likely benign", "likely pathogenic", "pathogenic")) %>% 
  ggplot(
    mapping = aes(
      x = factor(class, levels = c("benign", "likely benign", "likely pathogenic", "pathogenic")), 
      y = cosmis, 
      fill = factor(class, levels = c("benign", "likely benign", "likely pathogenic", "pathogenic"))
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
    values = pnw_palette(name = "Shuksan2", n = 7)[c(2, 3, 5, 6)]
  ) + 
  geom_abline(
    slope = 0, intercept = 0, color = "grey", size = 1, linetype = "dashed"
  ) +
  labs(
    y = "COSMIS score"
  ) +
  scale_x_discrete(
    labels = c("Benign","Likely Benign", "Likely Pathogenic", "Pathogenic")
  ) +
  scale_y_continuous(
    breaks = seq(-6, 6, 2),
    limits = c(-6, 6, 2),
    labels = sprintf("%.1f", seq(-6, 6, 2)),
    expand = c(0, 0)
  ) +
  plot_theme + theme(
    axis.title.x = element_blank()
  )

# save the violin plot to disk
ggsave(
  filename = paste(
    figure_dir, 
    "fig_sx_clinvar_cosmis_scores.svg", 
    sep = "/"
  ),
  plot = cosmis_violin_clinvar,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)