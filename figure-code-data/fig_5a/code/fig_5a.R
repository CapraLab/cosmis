rm(list = ls())
library(tidyverse)
library(PNWColors)

# data and figure paths
data_dir <- "/path/to/source_data/"

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

# load variant IDs
clinvar_patho_ids <-  read_tsv(
  file = paste(
    data_dir, "fig_5a/data",
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
    data_dir, "fig_5a/data",
    "clinvar_unambiguous_all_benign_ids_20210807.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)
clinvar_benign_ids <- clinvar_benign_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)
clinvar_vus_ids <-  read_tsv(
  file = paste(
    data_dir, "fig_5a/data",
    "clinvar_unambiguous_vus_ids_20210807.tsv", 
    sep = "/"
  ), 
  col_names = TRUE
)
clinvar_vus_ids <- clinvar_vus_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)

# cosmis scores for ClinVar pathogenic variants
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
vus_cosmis <- combined %>% filter(
  full_id %in% clinvar_vus_ids$full_id
) %>% mutate(
  label = 0.5
)

cosmis_clinvar <- rbind(
  patho_cosmis,
  benign_cosmis,
  vus_cosmis
)
cosmis_clinvar <- cosmis_clinvar %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", ifelse(label == 1, "pathogenic", "vus")), 
      levels = c("benign", "pathogenic", "vus")
    )
  )

# set plot theme
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
fig_5a <- cosmis_clinvar %>% 
  filter(class %in% c("benign", "pathogenic", "vus")) %>% 
  ggplot(
    mapping = aes(
      x = factor(class, levels = c("benign", "vus", "pathogenic")), 
      y = cosmis, 
      fill = factor(class, levels = c("benign", "vus", "pathogenic"))
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
    values = c(
      pnw_palette(name = "Shuksan2", n = 7)[2],
      "grey",
      pnw_palette(name = "Shuksan2", n = 7)[6]
    )
  ) + 
  geom_abline(
    slope = 0, intercept = 0, color = "grey", size = 1, linetype = "dashed"
  ) +
  labs(
    y = "COSMIS score"
  ) +
  scale_x_discrete(
    labels = c("Benign","VUS", "Pathogenic")
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
    data_dir, 
    "fig_5a/fig_5a.svg", 
    sep = "/"
  ),
  plot = fig_5a,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)