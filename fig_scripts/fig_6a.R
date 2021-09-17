# load required libraries
library(plotROC)
library(tidyverse)
library(PNWColors)
rm(list = ls())


# data and figure paths
data_dir <- "/path/to/datasets"
figure_dir <- "/path/to/where/figures/are/stored"

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
  full_id = str_c(uniprot_id, uniprot_pos, sep = "")
)

# load ClinVar variant IDs
clinvar_patho_ids <-  read_tsv(
  file = paste(
    result_dir, 
    "clinvar_unambiguous_pathogenic_ids_20210807.tsv", 
    sep = "/"
  ), 
  col_names = TRUE
)
clinvar_patho_ids <- clinvar_patho_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "")
)
clinvar_benign_ids <- read_tsv(
  file = paste(
    result_dir, 
    "clinvar_unambiguous_benign_ids_20210807.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)
clinvar_benign_ids <- clinvar_benign_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "")
)

patho_cosmis <- combined %>% filter(
  full_id %in% clinvar_patho_ids$full_id
) %>% mutate(
  label = 1
)
benign_cosmis <- combined %>% filter(
  full_id %in% clinvar_benign_ids$full_id
) %>% mutate(
  label = 0
)

cosmis_clinvar <- rbind(
  patho_cosmis,
  benign_cosmis
)
cosmis_clinvar <- cosmis_clinvar %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load pLI and Missense_Z scores
patho_pli_mis_z <- read_tsv(
  file = paste(
    data_dir,
    "clinvar_unambiguous_pathogenic_pli_mis_z_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_pli_mis_z$label <- 1
benign_pli_mis_z <- read_tsv(
  file = paste(
    data_dir,
    "clinvar_unambiguous_benign_pli_mis_z_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_pli_mis_z$label <- 0
pli_mis_z_clinvar <- bind_rows(
  patho_pli_mis_z,
  benign_pli_mis_z
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load MTR scores
patho_mtr3d <- read_tsv(
  file = paste(
    data_dir,
    "clinvar_unambiguous_pathogenic_mtr3d_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_mtr3d$label <- 1
benign_mtr3d <- read_tsv(
  file = paste(
    data_dir,
    "clinvar_unambiguous_benign_mtr3d_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_mtr3d$label <- 0
mtr3d_clinvar <- bind_rows(
  patho_mtr3d,
  benign_mtr3d
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

patho_mtr <- read_tsv(
  file = paste(
    data_dir,
    "clinvar_unambiguous_pathogenic_mtr_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_mtr$label <- 1
benign_mtr <- read_tsv(
  file = paste(
    data_dir,
    "clinvar_unambiguous_benign_mtr_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_mtr$label <- 0
mtr_clinvar <- bind_rows(
  patho_mtr,
  benign_mtr
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load RVIS scores
patho_rvis <- read_tsv(
  file = paste(
    data_dir,
    "clinvar_unambiguous_pathogenic_rvis_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_rvis$label <- 1
benign_rvis <- read_tsv(
  file = paste(
    data_dir,
    "clinvar_unambiguous_benign_rvis_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_rvis$label <- 0
rvis_clinvar <- bind_rows(
  patho_rvis,
  benign_rvis
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# COSMIS vs other metrics
cosmis_vs_others <- bind_cols(
  cosmis_clinvar %>% dplyr::select(uniprot_id, enst_id, uniprot_pos, uniprot_aa, cosmis), 
  mtr_clinvar %>% dplyr::select(mtr), mtr3d_clinvar %>% dplyr::select(mtr3d),
  rvis_clinvar %>% dplyr::select(rvis), pli_mis_z_clinvar %>% dplyr::select(pli, mis_z),
  cosmis_clinvar %>% dplyr::select(class)
) %>% drop_na()

# ggplot theme for ROC plot
roc_plot_theme <- theme_classic() + theme(
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

clinvar_cosmis_roc <- ggplot() +
  geom_roc(
    data = cosmis_vs_others,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = - cosmis
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 11)[11]
  ) +
  geom_roc(
    data = cosmis_vs_others,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = - mtr3d
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 11)[9]
  ) +
  geom_roc(
    data = cosmis_vs_others,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = - mtr
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 11)[7]
  ) +
  geom_roc(
    data = cosmis_vs_others,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = - rvis
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 11)[5]
  ) +
  geom_roc(
    data = cosmis_vs_others,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = pli
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 11)[3]
  ) +
  geom_roc(
    data = cosmis_vs_others,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = mis_z
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 11)[1]
  ) +
geom_abline(
  slope = 1, 
  intercept = 0, 
  linetype = "dashed", 
  size = 1
) + 
  xlab(
    "False positive rate"
  ) + 
  ylab(
    "True positive rate"
  ) + 
  scale_x_continuous(
    breaks = seq(0, 1, 0.2), 
    labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), 
    limits = c(0.0, 1.0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2), 
    labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), 
    limits = c(0.0, 1.0)
  ) +
  roc_plot_theme

# compute AUC values
# clinvar_cosmis_roc <- clinvar_cosmis_roc + 
#   annotate(
#     "text", 
#     x = 0.75, 
#     y = 0.25, 
#     size = 8,
#     label = paste("AUC =", round(calc_auc(clinvar_cosmis_roc)$AUC, 3))
#   )

# save the ROC plot to disk
ggsave(
  filename = paste(figure_dir, "cosmis_vs_other_metrics_roc.svg", sep = "/"),
  plot = clinvar_cosmis_roc,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)
