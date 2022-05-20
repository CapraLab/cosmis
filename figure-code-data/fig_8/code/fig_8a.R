# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())


# path to source data folder
data_dir <- "/path/to/source_data"

headers <- cols(
  uniprot_id = col_character(),
  enst_id = col_character(),
  uniprot_pos = col_integer(),
  uniprot_aa = col_character(),
  seq_separations = col_character(),
  num_contacts = col_integer(),
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

# load ion channel COSMIS scores
kcn_monomer_cosmis <- read_tsv(
  file = paste(data_dir, "fig_8/data", "KCN_monomer_cosmis.tsv", sep = "/"), 
  col_names = TRUE,
  col_types = headers
) %>% mutate(
  cosmis = (cs_mis_obs - mis_pmt_mean) / mis_pmt_sd,
  keys = str_c(uniprot_id, uniprot_pos, sep = "_")
)

kcn_multimer_cosmis <- read_tsv(
  file = paste(data_dir, "fig_8/data", "KCN_multimer_cosmis.tsv", sep = "/"), 
  col_names = TRUE,
  col_types = headers
) %>% mutate(
  cosmis = (cs_mis_obs - mis_pmt_mean) / mis_pmt_sd,
  keys = str_c(uniprot_id, uniprot_pos, sep = "_")
)

# load pathogenic and benign variants
pathogenic_ids <- read_tsv(
  file = paste(
    data_dir, "fig_8/data",
    "clinvar_KCN_unambiguous_pathogenic_uniprot_pos_20210807.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)
benign_ids <- read_tsv(
  file = paste(
    data_dir, "fig_8/data",
    "clinvar_KCN_unambiguous_benign_uniprot_pos_20210807.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)

# extract monomer and multimer COSMIS scores
interface_cosmis_mononer <- kcn_monomer_cosmis %>% filter(
  kcn_multimer_cosmis$num_contacts - kcn_monomer_cosmis$num_contacts > 0
) %>% mutate(
  class = "Interface"
)
noninterface_cosmis_mononer <- kcn_monomer_cosmis %>% filter(
  kcn_multimer_cosmis$num_contacts - kcn_monomer_cosmis$num_contacts == 0
) %>% mutate(
  class = "Non-interface"
)
interface_cosmis_multimer <- kcn_multimer_cosmis %>% filter(
  kcn_multimer_cosmis$num_contacts - kcn_monomer_cosmis$num_contacts > 0
) %>% mutate(
  class = "Interface*",
)
noninterface_cosmis_multimer <- kcn_multimer_cosmis %>% filter(
  kcn_multimer_cosmis$num_contacts - kcn_monomer_cosmis$num_contacts == 0
) %>% mutate(
  class = "Non-interface*"
)

# set up plot parameters
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
fig_8a <- bind_rows(
  interface_cosmis_mononer, noninterface_cosmis_mononer
) %>% 
  ggplot(
    mapping = aes(
      x = factor(class, levels = c("Non-interface", "Interface")), 
      y = cosmis, 
      fill = factor(class, levels = c("Non-interface", "Interface"))
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
    values = pnw_palette(name = "Shuksan2", n = 7)[c(2, 6)]
  ) + 
  labs(
    y = "COSMIS score"
  ) +
  scale_x_discrete(
    labels = c("Non-interface", "Interface")
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

ggsave(
  filename = paste(
    data_dir, "fig_8",
    "fig_8a.svg", 
    sep = "/"
  ),
  plot = fig_8a,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)
