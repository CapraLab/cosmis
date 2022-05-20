# load required packages
rm(list = ls())
library(tidyverse)
library(PNWColors)


# path to source data folder
data_dir <- "/path/to/source_data"

# load data set
prob_vs_count <- read_tsv(
  file = paste(
    data_dir, "fig_3a/data",
    "fig_3a_source_data_1.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)

# regress count of synonymous variants on synonymous mutability
syn_fit <- lm(
  formula = syn_count ~  syn_prob,
  data = prob_vs_count
)
prob_vs_count <- prob_vs_count %>% mutate(
  syn_exp = predict(
    object = syn_fit, 
    newdata = data.frame(syn_prob = syn_prob)
  ),
  mis_exp = predict(
    object = syn_fit,
    newdata = data.frame(syn_prob = mis_prob)
  )
)

# write data to disk
write_tsv(
  x = prob_vs_count,
  path = paste(
    data_dir, "fig_3a/data",
    "gnomad_all_enst_variant_counts.tsv",
    sep = "/"
  )
)

# Ensembl transcripts IDs for proteins with COSMIS scores
enst_with_cosmis <- read_tsv(
  file = paste(
    data_dir, 
    "fig_3a/data", 
    "fig_3a_source_data_2.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)
cosmis_prob_vs_count <- prob_vs_count %>% filter(
  enst_id %in% enst_with_cosmis$enst_id
)

# wilcoxon test
syn_dev <- cosmis_prob_vs_count$syn_count - cosmis_prob_vs_count$syn_exp
mis_dev <- cosmis_prob_vs_count$mis_count - cosmis_prob_vs_count$mis_exp
w_test <- wilcox.test(mis_dev, syn_dev, alternative = "less")

# make density plot
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
  aspect.ratio = 1.0,
  plot.margin = plot_margin
)

# make a density plot
fig_3a <- cosmis_prob_vs_count %>% ggplot() +
  geom_density(
    mapping = aes(x = syn_count - syn_exp),
    fill = pnw_palette(name = "Shuksan2", n = 5)[1],
    color = pnw_palette(name = "Shuksan2", n = 5)[1],
    alpha = 0.5,
    size = 1.2
  ) +
  geom_density(
    mapping = aes(
      x = mis_count - mis_exp
    ),
    fill = pnw_palette(name = "Shuksan2", n = 5)[5],
    color = pnw_palette(name = "Shuksan2", n = 5)[5],
    alpha = 0.5,
    size = 1.2
  ) +
  scale_x_continuous(
    name = "Observed - Expected",
    limits = c(-500, 500),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Density",
    limits = c(0, 0.030),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_3a", "fig_3a.svg", sep = "/"),
  plot = fig_3a,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)
