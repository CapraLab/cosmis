# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())

# path to source data folder
data_dir <- "/path/to/source_data"

# load data sets
cosmis_pdb <- read_tsv(
  file = paste(data_dir, "cosmis_scores_pdb.tsv.gz", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    cosmis = col_double(),
    cosmis_pvalue = col_double(),
    cossyn = col_double(),
    str_source = col_character()
  )
)
cosmis_swiss_model <- read_tsv(
  file = paste(data_dir, "cosmis_scores_swiss_model.tsv.gz", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    cosmis = col_double(),
    cosmis_pvalue = col_double(),
    cossyn = col_double(),
    str_source = col_character()
  )
)
cosmis_alphafold <- read_tsv(
  file = paste(data_dir, "cosmis_scores_alphafold.tsv.gz", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    cosmis = col_double(),
    cosmis_pvalue = col_double(),
    cossyn = col_double(),
    str_source = col_character()
  )
)
combined <- bind_rows(
  cosmis_pdb,
  cosmis_swiss_model,
  cosmis_alphafold
)
combined <- combined %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)


#===============================================================================
# load other scores
#===============================================================================
# gerp score
pdb_gerp <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_pdb_gerp_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_gerp <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_swiss_model_gerp_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_gerp <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_alphafold_gerp_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
gerp_scores <- bind_rows(
  pdb_gerp, swiss_model_gerp, alphafold_gerp
)

# phyloP scores
pdb_phylop <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_pdb_phylop_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_phylop <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_swiss_model_phylop_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_phylop <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_alphafold_phylop_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
phylop_scores <- bind_rows(
  pdb_phylop, swiss_model_phylop, alphafold_phylop
)

# phastCons scores
pdb_phastconst <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_pdb_phastcons_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_phastconst <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_swiss_model_phastcons_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_phastconst <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_alphafold_phastcons_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
phastcons_scores <- bind_rows(
  pdb_phastconst, swiss_model_phastconst, alphafold_phastconst
)

# ConSurf scores
pdb_consurf <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_pdb_consurf_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_consurf <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_swiss_model_consurf_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_consurf <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_alphafold_consurf_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
consurf_scores <- bind_rows(
  pdb_consurf, swiss_model_consurf, alphafold_consurf
)

# pLI missense_z
pdb_pli <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_pdb_pli_mis_z_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_pli <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_swiss_model_pli_mis_z_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_pli <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_alphafold_pli_mis_z_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
pli_scores <- bind_rows(
  pdb_pli, swiss_model_pli, alphafold_pli
)

# RVIS scores
pdb_rvis <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_pdb_rvis_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_rvis <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_swiss_model_rvis_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_rvis <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_alphafold_rvis_scores.tsv.gz", sep = "/"),
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
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_pdb_mtr_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_mtr <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_swiss_model_mtr_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_mtr <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_alphafold_mtr_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
mtr_scores <- bind_rows(
  pdb_mtr, swiss_model_mtr, alphafold_mtr
)

# MTR3D scores
pdb_mtr3d <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_pdb_mtr3d_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_mtr3d <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_swiss_model_mtr3d_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_mtr3d <- read_tsv(
  file = str_c(data_dir, "fig_6b/data", "cosmis_dataset_alphafold_mtr3d_scores.tsv.gz", sep = "/"),
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
  "gerp", "phylop", "phastcons", "consurf", "cosmis", "mtr3d", 
  "mtr", "rvis", "pli", "mis_z"
)
cors <- abs(cor(combined_na_removed[, scores], method = "spearman"))

# make a heatmap
fig_6b <- pheatmap(
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
  labels_row = c("GERP++", "phyloP", "phastCons", "R4S", "COSMIS",
                 "MTR", "RVIS", "pLI", "Missense_z"),
  labels_col = c("GERP++", "phyloP", "phastCons", "R4S", "COSMIS",
                 "MTR", "RVIS","pLI", "Misense_z"),
  fontsize = 16
)

ggsave(
  filename = paste(data_dir, "fig_6b", "fig_6b.svg", sep = "/"),
  plot = fig_6b,
  width = 5,
  height = 5,
  units = "in",
  device = "svg"
)