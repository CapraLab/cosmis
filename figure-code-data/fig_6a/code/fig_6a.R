# load required libraries
library(plotROC)
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

# load ClinVar variant IDs
clinvar_patho_ids <-  read_tsv(
  file = paste(
    data_dir, "fig_6a/data",
    "clinvar_unambiguous_all_pathogenic_ids_20210807.tsv", 
    sep = "/"
  ), 
  col_names = TRUE
)
clinvar_patho_ids <- clinvar_patho_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)
clinvar_benign_ids <- read_tsv(
  file = paste(
    data_dir, "fig_6a/data",
    "clinvar_unambiguous_all_benign_ids_20210807.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)
clinvar_benign_ids <- clinvar_benign_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
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
    data_dir, "fig_6a/data",
    "clinvar_unambiguous_pathogenic_pli_mis_z_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_pli_mis_z$label <- 1
benign_pli_mis_z <- read_tsv(
  file = paste(
    data_dir, "fig_6a/data",
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
    data_dir, "fig_6a/data",
    "clinvar_unambiguous_pathogenic_mtr3d_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_mtr3d$label <- 1
benign_mtr3d <- read_tsv(
  file = paste(
    data_dir, "fig_6a/data",
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
    data_dir, "fig_6a/data",
    "clinvar_unambiguous_pathogenic_mtr_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_mtr$label <- 1
benign_mtr <- read_tsv(
  file = paste(
    data_dir, "fig_6a/data",
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
    data_dir, "fig_6a/data",
    "clinvar_unambiguous_pathogenic_rvis_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_rvis$label <- 1
benign_rvis <- read_tsv(
  file = paste(
    data_dir, "fig_6a/data",
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

fig_6a <- ggplot() +
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
  filename = paste(data_dir, "fig_6a", "fig_6a.svg", sep = "/"),
  plot = fig_6a,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)