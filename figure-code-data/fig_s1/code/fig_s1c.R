# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())

# path to source data folder
data_dir <- "/path/to/source_data"

# load alphafold dataset
cosmis_df_alphafold <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_alphafold.tsv.gz", sep = "/"),
  col_names = TRUE,
  col_types = cols(
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
    enst_mis_exp = col_double(),
    plddt = col_double(),
    uniprot_length = col_integer()
  )
)

# compute median pLDDT scores
median_plddt <- cosmis_df_alphafold %>% group_by(uniprot_id) %>% summarise(
  median_plddt = median(plddt)
)

# make a violin plot
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme_classic() + theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text = element_text(
    size = 20, color = "black"
  ),
  axis.text.x = element_text(
    margin = margin(t = 10)
  ),
  axis.title.x = element_text(
    color = "black", size = 24, margin = margin(t = 10)
  ),
  axis.title.y = element_text(
    color = "black", size = 24, margin = margin(r = 10)
  ),
  axis.ticks.length.y = unit(.25, "cm"),
  axis.ticks.y = element_line(),
  axis.ticks.x = element_blank(),
  legend.position = "none",
  plot.margin = plot_margin
)

fig_s1c <- median_plddt %>% 
  ggplot(
    mapping = aes(
      x = "", 
      y = median_plddt
    )
  ) +
  geom_violin(
    trim = FALSE,
    size = 1,
    fill = pnw_palette(name = "Shuksan2", n = 7)[c(6)]
  ) +
  geom_boxplot(
    width = 0.1,
    fill = "white",
    outlier.shape = NA
  ) +
  labs(
    x = "",
    y = str_wrap("Median pLDDT of AF2 models", width = 30)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 100), 
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s1", "fig_s1c.svg", sep = "/"),
  plot = fig_s1c,
  width = 4,
  height = 6,
  units = "in",
  device = "svg",
)
