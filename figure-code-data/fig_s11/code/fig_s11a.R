library(tidyverse)
library(PNWColors)
library(plotROC)
rm(list = ls())

# load data sets
data_dir <-  "/Users/lib14/code_space/cosmis/source_data"

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

scn_cosmis <- read_tsv(
  file = paste(data_dir, "fig_s11/data", "SCN_cosmis.tsv", sep = "/"), 
  col_names = TRUE,
  col_types = headers
) %>% mutate(
  cosmis = (cs_mis_obs - mis_pmt_mean) / mis_pmt_sd,
  keys = str_c(uniprot_id, uniprot_pos, sep = "_")
)
scn_mtr <- read_csv(
  file = paste(data_dir, "fig_s11/data", "SCN_mtr_by_uniprot.csv", sep = "/"),
  col_names = TRUE
)

# SCN5A functional data
glazer_benign_ids <- read.table(
  file = paste(
    data_dir,
    "fig_s11/data",
    "glazer_benign_ids.txt", 
    sep = "/"
  ),
  header = TRUE
)
glazer_pathogenic_ids <- read.table(
  file = paste(
    data_dir,
    "fig_s11/data",
    "glazer_pathogenic_ids.txt", 
    sep = "/"
  ),
  header = TRUE
)

# create datasets
cosmis_pathogenic <- scn_cosmis %>% filter(
  keys %in% glazer_pathogenic_ids$pos
) %>% mutate(
  label = 1,
  class = "pathogenic"
)
cosmis_benign <- scn_cosmis %>% filter(
  keys %in% glazer_benign_ids$pos
) %>% mutate(
  label = 0,
  class = "benign"
)
glazer_cosmis <- bind_rows(
  cosmis_benign, cosmis_pathogenic
)

# create datasets
mtr_pathogenic <- scn_mtr %>% filter(
  uniprot_aapos %in% glazer_pathogenic_ids$pos
) %>% mutate(
  label = 1,
  class = "pathogenic"
)
mtr_benign <- scn_mtr %>% filter(
  uniprot_aapos %in% glazer_benign_ids$pos
) %>% mutate(
  label = 0,
  class = "benign"
)
glazer_mtr <- bind_rows(
  mtr_benign, mtr_pathogenic
)

# ggplot theme for ROC plot
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
roc_plot_theme <- theme_classic() + theme(
  panel.border = element_rect(
    colour = "black", 
    size = 1, 
    fill = "transparent"
  ),
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
  aspect.ratio = 1.0,
  plot.margin = plot_margin
)

fig_s11a <- ggplot() +
  geom_roc(
    data = glazer_cosmis,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = - cosmis
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 7)[7]
  ) +
  geom_roc(
    data = glazer_mtr,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = - mtr 
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 7)[5]
  ) +
  geom_abline(
    slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    size = 1
  ) + 
  xlab(
    "False positive rate"
  ) + 
  ylab(
    "True positive rate"
  ) + 
  scale_x_continuous(
    breaks = seq(0, 1, 0.2), 
    labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), 
    limits = c(0.0, 1.0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2), 
    labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), 
    limits = c(0.0, 1.0)
  ) +
  roc_plot_theme

ggsave(
  filename = paste(data_dir, "fig_s11", "fig_s11a.svg", sep = "/"),
  plot = fig_s11a,
  width = 6,
  height = 6,
  units = "in",
  device = "svg",
)