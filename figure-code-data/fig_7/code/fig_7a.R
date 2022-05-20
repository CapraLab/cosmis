# load required libraries
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
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_pdb_gerp_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_gerp <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_swiss_model_gerp_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_gerp <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_alphafold_gerp_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
gerp_scores <- bind_rows(
  pdb_gerp, swiss_model_gerp, alphafold_gerp
)

# phyloP scores
pdb_phylop <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_pdb_phylop_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_phylop <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_swiss_model_phylop_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_phylop <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_alphafold_phylop_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
phylop_scores <- bind_rows(
  pdb_phylop, swiss_model_phylop, alphafold_phylop
)

# phastCons scores
pdb_phastconst <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_pdb_phastcons_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_phastconst <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_swiss_model_phastcons_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_phastconst <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_alphafold_phastcons_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
phastcons_scores <- bind_rows(
  pdb_phastconst, swiss_model_phastconst, alphafold_phastconst
)

# ConSurf scores
pdb_consurf <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_pdb_consurf_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_consurf <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_swiss_model_consurf_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_consurf <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_alphafold_consurf_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
consurf_scores <- bind_rows(
  pdb_consurf, swiss_model_consurf, alphafold_consurf
)

# pLI missense_z
pdb_pli <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_pdb_pli_mis_z_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_pli <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_swiss_model_pli_mis_z_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_pli <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_alphafold_pli_mis_z_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
pli_scores <- bind_rows(
  pdb_pli, swiss_model_pli, alphafold_pli
)

# RVIS scores
pdb_rvis <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_pdb_rvis_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_rvis <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_swiss_model_rvis_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_rvis <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_alphafold_rvis_scores.tsv.gz", sep = "/"),
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
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_pdb_mtr_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_mtr <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_swiss_model_mtr_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_mtr <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_alphafold_mtr_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
mtr_scores <- bind_rows(
  pdb_mtr, swiss_model_mtr, alphafold_mtr
)

# MTR3D scores
pdb_mtr3d <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_pdb_mtr3d_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
swiss_model_mtr3d <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_swiss_model_mtr3d_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
alphafold_mtr3d <- read_tsv(
  file = str_c(data_dir, "fig_7/data", "cosmis_dataset_alphafold_mtr3d_scores.tsv.gz", sep = "/"),
  col_names = TRUE
)
mtr3d_scores <- bind_rows(
  pdb_mtr3d, swiss_model_mtr3d, alphafold_mtr3d
)


combined_na_removed <- bind_cols(
  combined %>% dplyr::select(uniprot_id, enst_id, uniprot_pos, uniprot_aa, cosmis), 
  gerp_scores %>% dplyr::select(gerp), phylop_scores %>% dplyr::select(phylop), 
  phastcons_scores %>% dplyr::select(phastcons), consurf_scores %>% dplyr::select(consurf), 
  pli_scores %>% dplyr::select(pli, mis_z),
  mtr_scores %>% dplyr::select(mtr), rvis_scores %>% dplyr::select(rvis), 
  mtr3d_scores %>% dplyr::select(mtr3d)
) %>% drop_na()

# load de novo variants coordinates
case_coords <- read_table(
  file = str_c(
    data_dir, "fig_7/data", "de_novo_case_variant_coord.txt", sep = "/"
  ),
  col_names = TRUE
)
control_coords <- read_table(
  file = str_c(
    data_dir, "fig_7/data", "de_novo_control_variant_coord.txt", sep = "/"
  ),
  col_names = TRUE
)

combined_na_removed <- combined_na_removed %>% mutate(
  keys = str_c(enst_id, uniprot_pos, sep = "_")
)

control_cosmis_df <- combined_na_removed %>% filter(
  keys %in% control_coords$coord
)
case_cosmis_df <- combined_na_removed %>% filter(
  keys %in% case_coords$coord
)

case_cosmis_df <- case_cosmis_df %>% mutate(
  class = "case"
)
control_cosmis_df <- control_cosmis_df %>% mutate(
  class = "control"
)

de_novo_cosmis <- bind_rows(
  case_cosmis_df, control_cosmis_df
) %>% dplyr::select(-keys)

# make a violin plot
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
  plot.margin = plot_margin
)

# make a violin plot
fig_7a <- ggplot(
  data = de_novo_cosmis,
  mapping = aes(
    x = factor(class, levels = c("control", "case")), 
    y = cosmis, 
    fill = factor(class, levels = c("control", "case"))
  )
) +
  geom_violin(
    trim = FALSE,
    size = 1,
    alpha = 0.5
  ) +
  geom_boxplot(
    width = 0.1,
    fill = "white",
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = pnw_palette(name = "Shuksan2", n = 7)[c(2, 6)]
  ) + 
  labs(
    y = "COSMIS score"
  ) +
  scale_x_discrete(
    labels = c("Control", "Case")
  ) +
  scale_y_continuous(
    breaks = seq(-6, 6, 2),
    limits = c(-6, 6), 
    labels = sprintf("%.1f", seq(-6, 6, 2)),
    expand = c(0, 0)
  ) +
  plot_theme + theme(
    axis.title.x = element_blank()
  )

ggsave(
  filename = paste(
    data_dir, "fig_7", 
    "fig_7a.svg", 
    sep = "/"
  ),
  plot = fig_7a,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)
