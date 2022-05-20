# load required packages
library(tidyverse)
library(PNWColors)
library(plotROC)
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

# create datasets
pathogenic_multimer_df <- kcn_multimer_cosmis %>% filter(
  keys %in% pathogenic_ids$id
) %>% mutate(
  label = 1,
  class = "pathogenic"
)
benign_multimer_df <- kcn_multimer_cosmis %>% filter(
  keys %in% benign_ids$id
) %>% mutate(
  label = 0,
  class = "benign"
)

# plot theme
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

fig_8c <- ggplot() +
  geom_roc(
    data = bind_rows(pathogenic_multimer_df, benign_multimer_df),
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = - cosmis
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.5,
    colour = pnw_palette(name = "Shuksan2", n = 3)[3]
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

# compute AUC values
kcn_cosmis_roc <- kcn_cosmis_roc + 
  annotate(
    "text", 
    x = 0.75, 
    y = 0.25, 
    size = 8,
    label = paste("AUC =", round(calc_auc(kcn_cosmis_roc)$AUC, 3))
  )

# save the ROC plot to disk
ggsave(
  filename = paste(data_dir, "fig_8", "fig_8c.svg", sep = "/"),
  plot = fig_8c,
  width = 6,
  height = 6,
  units = "in",
  device = "svg",
)
