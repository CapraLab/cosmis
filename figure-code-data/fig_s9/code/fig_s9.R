library(tidyverse)
library(PNWColors)
rm(list = ls())

# path to source data folder
data_dir <- "/path/to/source_data"

## load data sets
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
combined <- rbind(
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
    data_dir, "fig_s9/data",
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
    data_dir, "fig_s9/data",
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
    data_dir, "fig_s9/data",
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
    data_dir, "fig_s9/data",
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
fig_s9 <- cosmis_clinvar %>% 
  ggplot(
    mapping = aes(
      x = class,
      y = cosmis, 
      fill = class
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
    data_dir, "fig_s9", 
    "fig_s9.svg", 
    sep = "/"
  ),
  plot = fig_s9,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)