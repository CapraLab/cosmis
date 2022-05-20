# load required packages
rm(list = ls())
library(tidyverse)
library(PNWColors)
library(plotROC)


# load data sets
data_dir <-  "/path/to/source_data"

headers <-  cols(
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
  uniprot_length = col_integer(),
  cosmis = col_double()
)


# load INSIDER COSMIS scores
insider_monomer_cosmis <- read_tsv(
  file = paste(
    data_dir, "fig_s14/data", "cosmis_dataset_insider_monomer.tsv.gz", sep = "/"
  ), 
  col_names = TRUE,
  col_types = headers
) %>% mutate(
  keys = str_c(uniprot_id, uniprot_pos, sep = "_")
)

insider_multimer_cosmis <- read_tsv(
  file = paste(
    data_dir, "fig_s14/data", "cosmis_dataset_insider_multimer.tsv.gz", sep = "/"
  ), 
  col_names = TRUE,
  col_types = headers
) %>% mutate(
  keys = str_c(uniprot_id, uniprot_pos, sep = "_")
)

# extract monomer and multimer COSMIS scores
insider_multimer_cosmis <- insider_multimer_cosmis %>% filter(
  insider_monomer_cosmis$num_contacts < 11
)
insider_monomer_cosmis <- insider_monomer_cosmis %>% filter(
  num_contacts < 11
)

interface_cosmis_mononer <- insider_monomer_cosmis %>% filter(
  insider_multimer_cosmis$num_contacts - insider_monomer_cosmis$num_contacts > 0
) %>% mutate(
  class = "Interface"
)
noninterface_cosmis_mononer <- insider_monomer_cosmis %>% filter(
  insider_multimer_cosmis$num_contacts - insider_monomer_cosmis$num_contacts == 0
) %>% mutate(
  class = "Non-interface"
)
interface_cosmis_multimer <- insider_multimer_cosmis %>% filter(
  insider_multimer_cosmis$num_contacts - insider_monomer_cosmis$num_contacts > 0
) %>% mutate(
  class = "Interface",
)
noninterface_cosmis_multimer <- insider_multimer_cosmis %>% filter(
  insider_multimer_cosmis$num_contacts - insider_monomer_cosmis$num_contacts == 0
) %>% mutate(
  class = "Non-interface"
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
fig_s14a <- bind_rows(
    interface_cosmis_multimer, noninterface_cosmis_multimer
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
  geom_abline(
    intercept = 0, slope = 0, color = "grey", linetype = "dashed"
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
    data_dir, "fig_s14", 
    "fig_s14a.svg", 
    sep = "/"
  ),
  plot = fig_s14a,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)

#===============================================================================
# End of multimer-based analysis
#===============================================================================

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

fig_s14b <-  ggplot() + 
  geom_histogram(
    mapping = aes(interface_cosmis_multimer$cosmis - interface_cosmis_mononer$cosmis),
    bins = 50,
    colour = "black",
    fill = pnw_palette(name = "Shuksan2", n = 7)[2]
  ) +
  geom_vline(
    xintercept = 0, color = "grey", linetype = "dashed"
  ) +
  scale_x_continuous(
    name = str_wrap(
      "COSMIS score difference (homodimer - monomer)", width = 30
    ),
    limits = c(-2.5, 2.5),
    breaks = seq(-2.5, 2.5, 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Count",
    limits = c(0, 5000),
    breaks = seq(0, 5000, 1000),
    expand = c(0, 0)
  ) +
  hist_plot_theme

ggsave(
  filename = paste(data_dir, "fig_s14", "fig_s14b.svg", sep = "/"),
  plot = fig_s14b,
  width = 8,
  height = 6,
  units = "in",
  device = "svg",
)