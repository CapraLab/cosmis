# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())

# path to source data folder
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

# load pLI scores
pli_pdb <- read_csv(
  file = paste(data_dir, "fig_4a/data/pdb_dataset_pli_mis_z.csv", sep = "/"),
  col_names = TRUE
)
pli_sm <- read_csv(
  file = paste(data_dir, "fig_4a/data/sm_dataset_pli_mis_z.csv", sep = "/"),
  col_names = TRUE
)
pli_af2 <- read_csv(
  file = paste(data_dir, "fig_4a/data/af2_dataset_pli_mis_z.csv", sep = "/"),
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
    mean_cosmis = mean(cosmis),
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

fig_4a <- cosmis_vs_pli %>% 
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
  filename = paste(data_dir, "fig_4a", "fig_4a.svg", sep = "/"),
  plot = fig_4a,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)