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

# create datasets
interface_cosmis_mononer <- kcn_monomer_cosmis %>% filter(
  kcn_multimer_cosmis$num_contacts - kcn_monomer_cosmis$num_contacts > 0
) %>% mutate(
  class = "Interface"
)
interface_cosmis_multimer <- kcn_multimer_cosmis %>% filter(
  kcn_multimer_cosmis$num_contacts - kcn_monomer_cosmis$num_contacts > 0
) %>% mutate(
  class = "Interface*",
)

# set plot parameters
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
hist_plot_theme <- theme_classic() + theme(
  panel.border = element_blank(),
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
  axis.ticks = element_line(),
  legend.position = "none",
  plot.margin = plot_margin
)

fig_8b <-  ggplot() + 
  geom_histogram(
    mapping = aes(interface_cosmis_multimer$cosmis - interface_cosmis_mononer$cosmis),
    bins = 40,
    colour = "black",
    fill = pnw_palette(name = "Shuksan2", n = 7)[2]
  ) +
  scale_x_continuous(
    name = str_wrap(
      "COSMIS score difference (oligomer - monomer)", width = 30
    ),
    limits = c(-2, 2),
    breaks = seq(-2, 2, 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Count",
    limits = c(0, 1000),
    breaks = seq(0, 1000, 200),
    expand = c(0, 0)
  ) +
  hist_plot_theme

ggsave(
  filename = paste(data_dir, "fig_8", "fig_8b.svg", sep = "/"),
  plot = fig_8b,
  width = 10,
  height = 6,
  units = "in",
  device = "svg",
)
