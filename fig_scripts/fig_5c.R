# load required packages
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
  cosmis = (cs_mis_obs - mis_pmt_mean) / mis_pmt_sd
)

# load haploinsufficient gene transcript IDs
hi_uniprot<- read_tsv(
  file = str_c(
    data_dir,
    "clingen_level3_haploinsufficient_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# load olfactory receptor gene transcript IDs
or_uniprot <- read_tsv(
  file = str_c(
    data_dir,
    "olfactory_receptor_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# load all dominant gene transcript IDs
dominant_uniprot <- read_tsv(
  file = str_c(
    data_dir,
    "all_dominant_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# load all dominant gene transcript IDs
recessive_uniprot <- read_tsv(
  file = str_c(
    data_dir,
    "all_recessive_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# load all dominant gene transcript IDs
essential_uniprot <- read_tsv(
  file = str_c(
    data_dir,
    "essential_genes_uniprot.tsv",
    sep = "/"
  ),
  col_names = TRUE
)

# load all dominant gene transcript IDs
non_essential_uniprot <- read_tsv(
  file = str_c(
    data_dir,
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

gene_class_cosmis_violin <- bind_rows(
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
  filename = paste(figure_dir, "gene_class_cosmis_violin_plot.svg", sep = "/"),
  plot = gene_class_cosmis_violin,
  width = 8,
  height = 10,
  units = "in",
  device = "svg",
)

