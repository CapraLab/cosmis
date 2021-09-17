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
combined <- combined %>% mutate(
  cosmis = (cs_mis_obs - mis_pmt_mean) / mis_pmt_sd,
  full_id = str_c(uniprot_id, uniprot_pos, sep = "")
)

#===============================================================================
# load other scores
#===============================================================================
# gerp score
pdb_gerp <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_gerp_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_gerp <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_gerp_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_gerp <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_gerp_scores.tsv", sep = "/"),
  col_names = TRUE
)
gerp_scores <- bind_rows(
  pdb_gerp, swiss_model_gerp, alphafold_gerp
)

# phyloP scores
pdb_phylop <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_phylop_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_phylop <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_phylop_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_phylop <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_phylop_scores.tsv", sep = "/"),
  col_names = TRUE
)
phylop_scores <- bind_rows(
  pdb_phylop, swiss_model_phylop, alphafold_phylop
)

# phastCons scores
pdb_phastconst <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_phastcons_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_phastconst <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_phastcons_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_phastconst <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_phastcons_scores.tsv", sep = "/"),
  col_names = TRUE
)
phastcons_scores <- bind_rows(
  pdb_phastconst, swiss_model_phastconst, alphafold_phastconst
)

# ConSurf scores
pdb_consurf <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_consurf_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_consurf <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_consurf_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_consurf <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_consurf_scores.tsv", sep = "/"),
  col_names = TRUE
)
consurf_scores <- bind_rows(
  pdb_consurf, swiss_model_consurf, alphafold_consurf
)

# pLI missense_z
pdb_pli <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_pli_mis_z_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_pli <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_pli_mis_z_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_pli <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_pli_mis_z_scores.tsv", sep = "/"),
  col_names = TRUE
)
pli_scores <- bind_rows(
  pdb_pli, swiss_model_pli, alphafold_pli
)

# RVIS scores
pdb_rvis <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_rvis_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_rvis <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_rvis_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_rvis <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_rvis_scores.tsv", sep = "/"),
  col_names = TRUE,
  col_types = cols(
    hgnc = col_character(),
    rvis = col_double(),
    rvis_percentile = col_double()
  )
)
rvis_scores <- bind_rows(
  pdb_rvis, swiss_model_rvis, alphafold_rvis
)

# MTR scores
pdb_mtr <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_mtr_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_mtr <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_mtr_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_mtr <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_mtr_scores.tsv", sep = "/"),
  col_names = TRUE
)
mtr_scores <- bind_rows(
  pdb_mtr, swiss_model_mtr, alphafold_mtr
)

# MTR3D scores
pdb_mtr3d <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_mtr3d_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_mtr3d <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_mtr3d_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_mtr3d <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_mtr3d_scores.tsv", sep = "/"),
  col_names = TRUE
)
mtr3d_scores <- bind_rows(
  pdb_mtr3d, swiss_model_mtr3d, alphafold_mtr3d
)

#===============================================================================
# correlations between different metrics
#===============================================================================
library(pheatmap)

combined_na_removed <- bind_cols(
  combined %>% dplyr::select(cosmis), gerp_scores %>% dplyr::select(gerp),
  phylop_scores %>% dplyr::select(phylop), phastcons_scores %>% dplyr::select(phastcons),
  consurf_scores %>% dplyr::select(consurf), pli_scores %>% dplyr::select(pli, mis_z),
  mtr_scores %>% dplyr::select(mtr), rvis_scores %>% dplyr::select(rvis), 
  mtr3d_scores %>% dplyr::select(mtr3d)
) %>% drop_na()

scores <- c(
  "gerp", "phylop", "phastcons", "consurf", "z_score", "mtr3d", 
  "mtr", "rvis", "pli", "mis_z"
)
cors <- abs(cor(combined_na_removed[, scores], method = "spearman"))

# make a heatmap
cor_heatmap <- pheatmap(
  cors, 
  border_color = "white", 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  labels_row = "",
  labels_col = "",
  display_numbers = TRUE,
  # legend = FALSE,
  legend_breaks = seq(0, 1.0, 0.2),
  legend_labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
  # labels_row = c("GERP++", "phyloP", "phastCons", "R4S", "COSMIS", 
  #                "MTR", "RVIS", "pLI", "Missense_z"),
  # labels_col = c("GERP++", "phyloP", "phastCons", "R4S", "COSMIS", 
  #                "MTR", "RVIS","pLI", "Misense_z"),
  fontsize = 16
)

ggsave(
  filename = paste(figure_dir, "constraint_score_correlation_heatmap.svg", sep = "/"),
  plot = cor_heatmap,
  width = 5,
  height = 5,
  units = "in",
  device = "svg"
)