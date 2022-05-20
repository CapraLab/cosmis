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

# load haploinsufficient gene transcript IDs
hi_uniprot<- read_tsv(
  file = str_c(
    data_dir, "fig_5c/data",
    "clingen_level3_haploinsufficient_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# load olfactory receptor gene transcript IDs
or_uniprot <- read_tsv(
  file = str_c(
    data_dir, "fig_5c/data",
    "olfactory_receptor_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# load all dominant gene transcript IDs
dominant_uniprot <- read_tsv(
  file = str_c(
    data_dir, "fig_5c/data",
    "all_dominant_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# load all dominant gene transcript IDs
recessive_uniprot <- read_tsv(
  file = str_c(
    data_dir, "fig_5c/data",
    "all_recessive_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# load all dominant gene transcript IDs
essential_uniprot <- read_tsv(
  file = str_c(
    data_dir, "fig_5c/data",
    "essential_genes_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# load all dominant gene transcript IDs
non_essential_uniprot <- read_tsv(
  file = str_c(
    data_dir, "fig_5c/data",
    "non_essential_genes_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# extract rows from COSMIS dataset
hi_cosmis <- combined %>% filter(
  uniprot_id %in% hi_uniprot$uniprot_id
) %>% mutate(
  class = "Haploinsufficient"
)

or_cosmis <- combined %>% filter(
  uniprot_id %in% or_uniprot$uniprot_id
) %>% mutate(
  class = "Olfactory"
)

dominant_cosmis <- combined %>% filter(
  uniprot_id %in% dominant_uniprot$uniprot_id
) %>% mutate(
  class = "Dominant"
)

recessive_cosmis <- combined %>% filter(
  uniprot_id %in% recessive_uniprot$uniprot_id
) %>% mutate(
  class = "Recessive"
)

essential_cosmis <- combined %>% filter(
  uniprot_id %in% essential_uniprot$uniprot_id
) %>% mutate(
  class = "Essential"
)

non_essential_cosmis <- combined %>% filter(
  uniprot_id %in% non_essential_uniprot$uniprot_id
) %>% mutate(
  class = "Non essential"
)

combined$class <- "All"

# make a violin plot
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.x = element_line(size = 1),
  axis.line.y = element_line(size = 1),
  axis.text = element_text(
    size = 24, color = "black"
  ),
  axis.title.x = element_text(
    color = "black", size = 28, margin = margin(t = 10)
  ),
  axis.title.y = element_text(
    color = "black", size = 28, margin = margin(r = 10)
  ),
  axis.ticks.length = unit(0.4, "cm"),
  axis.ticks = element_line(size = 1),
  legend.position = "none",
  plot.margin = plot_margin
)

fig_5c <- bind_rows(
  hi_cosmis, essential_cosmis, dominant_cosmis, combined, 
  recessive_cosmis, non_essential_cosmis, or_cosmis
) %>% mutate(
  class = factor(
    class, 
    levels = c(
      "Haploinsufficient", "Essential", "Dominant", "All", "Recessive", 
      "Non essential", "Olfactory"
    ),
    ordered = TRUE
  )
) %>% 
  ggplot(
    mapping = aes(
      x = class,
      y = cosmis,
      fill = class
    )
  ) +
  geom_hline(
    yintercept = 0, 
    linetype = "dotted", 
    color = "gray", 
    size = 1.5
  ) +
  geom_violin(
    size = 0.8
  ) +
  geom_boxplot(
    fill = "white",
    size = 0.8,
    width = 0.12,
    outlier.shape = NA
  ) +
  coord_flip(
  ) +
  scale_fill_manual(
    values = rev(pnw_palette(name = "Shuksan2", n = 7))
  ) +
  xlab(
    label = NULL
  ) +
  scale_y_continuous(
    name = str_wrap("COSMIS score", width = 30),
    limits = c(-6, 6),
    breaks = seq(-6, 6, 2),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(
    data_dir, "fig_5c", "fig_5c.svg", sep = "/"),
  plot = fig_5c,
  width = 8,
  height = 10,
  units = "in",
  device = "svg",
)