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

# fraction of constrained sites
by_uniprot <- combined %>% group_by(uniprot_id) %>% summarise(
  n_res = n(),
  num_cst_sites = mean(sum(cosmis_pvalue < 0.01)),
  frac_cst_sites = mean(sum(cosmis_pvalue < 0.01) / n_res),
)

# load haploinsufficient gene transcript IDs
hi_uniprot<- read_tsv(
  file = str_c(
    data_dir, "fig_s10/data",
    "clingen_level3_haploinsufficient_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

hi_df <- by_uniprot %>% filter(
  uniprot_id %in% hi_uniprot$uniprot_id
) %>% mutate(
  class = "Haploinsufficient"
)

# load olfactory receptor gene transcript IDs
or_uniprot <- read_tsv(
  file = str_c(
    data_dir, "fig_s10/data",
    "olfactory_receptor_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

or_df <- by_uniprot %>% filter(
  uniprot_id %in% or_uniprot$uniprot_id
) %>% mutate(
  class = "Olfactory"
)

# load all dominant gene transcript IDs
dominant_uniprot <- read_tsv(
  file = str_c(
    data_dir, "fig_s10/data",
    "all_dominant_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

dominant_df <- by_uniprot %>% filter(
  uniprot_id %in% dominant_uniprot$uniprot_id
) %>% mutate(
  class = "Dominant"
)

# load all dominant gene transcript IDs
recessive_uniprot <- read_tsv(
  file = str_c(
    data_dir, "fig_s10/data",
    "all_recessive_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

recessive_df <- by_uniprot %>% filter(
  uniprot_id %in% recessive_uniprot$uniprot_id
) %>% mutate(
  class = "Recessive"
)

# load all dominant gene transcript IDs
essential_uniprot <- read_tsv(
  file = str_c(
    data_dir, "fig_s10/data",
    "essential_genes_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

essential_df <- by_uniprot %>% filter(
  uniprot_id %in% essential_uniprot$uniprot_id
) %>% mutate(
  class = "Essential"
)

# load all dominant gene transcript IDs
non_essential_uniprot <- read_tsv(
  file = str_c(
    data_dir, "fig_s10/data",
    "non_essential_genes_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

non_essential_df <- by_uniprot %>% filter(
  uniprot_id %in% non_essential_uniprot$uniprot_id
) %>% mutate(
  class = "Non essential"
)

# combine gene categories
by_uniprot$class <- "All"
cst_fraction <- bind_rows(
  hi_df, or_df, dominant_df, recessive_df, essential_df, non_essential_df, by_uniprot
) %>% group_by(class) %>% summarise(
  fraction = sum(num_cst_sites > 0) / n()
) 

# make a bar plot
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text = element_text(
    size = 20, color = "black"
  ),
  axis.title.x = element_text(
    color = "black", size = 24, margin = margin(t = 10)
  ),
  axis.title.y = element_text(
    color = "black", size = 24, margin = margin(r = 10)
  ),
  axis.ticks.length = unit(.25, "cm"),
  axis.ticks = element_line(),
  legend.position = "none",
  plot.margin = plot_margin
)

fig_s10a <- ggplot(
    data = cst_fraction, 
    mapping = aes(
      x = reorder(class, fraction), 
      y = fraction,
      fill = reorder(class, fraction)
    )
  ) +
  geom_bar(
    stat = "identity"
  ) +
  coord_flip(
  ) +
  scale_fill_manual(
    values = pnw_palette(name = "Shuksan2", n = 7)
  ) + 
  xlab(
    label = NULL
  ) +
  scale_y_continuous(
    name = str_wrap("Fraction of genes with as least one high-confidence constraint site ", width = 35),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(
    data_dir, "fig_s10", "fig_s10a.svg", sep = "/"
  ),
  plot = fig_s10a,
  width = 8,
  height = 6,
  units = "in",
  device = "svg",
)

# make a boxplot
plot_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text = element_text(
    size = 20, color = "black"
  ),
  axis.title.x = element_text(
    color = "black", size = 24, margin = margin(t = 10)
  ),
  axis.title.y = element_text(
    color = "black", size = 24, margin = margin(r = 10)
  ),
  axis.ticks.length = unit(.25, "cm"),
  axis.ticks = element_line(),
  legend.position = "none",
  plot.margin = plot_margin
)

fig_s10b <- bind_rows(
  hi_df, essential_df, dominant_df, by_uniprot, recessive_df, non_essential_df, or_df
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
      y = frac_cst_sites,
      fill = class
    )
  ) +
  geom_boxplot(
    outlier.alpha = 0.1
  ) +
  scale_fill_manual(
    values = rev(pnw_palette(name = "Shuksan2", n = 7))
  ) +
  xlab(
    label = NULL
  ) +
  scale_y_continuous(
    name = str_wrap("Fraction of high-confidence constrained sites in each gene", width = 30),
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, 0.2),
    expand = c(0, 0)
  ) +
  plot_theme + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
    # axis.ticks.x = element_blank()
  )

ggsave(
  filename = paste(
    data_dir, "fig_s10", "fig_s10b.svg", sep = "/"
  ),
  plot = fig_s10b,
  width = 12,
  height = 7,
  units = "in",
  device = "svg",
)
