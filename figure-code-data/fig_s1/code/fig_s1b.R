# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())


# path to source data folder
data_dir <- "/path/to/source_data"

swiss_model_stat <- read_csv(
  file = paste(data_dir, "fig_s1/data", "fig_s1b_source_data.csv", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    iso_id = col_character(),
    uniprot_length = col_integer(),
    source = col_character(),
    from = col_integer(),
    to = col_integer(),
    template = col_character(),
    qmean = col_double(),
    qnorm = col_double(),
    seq_id = col_double(),
    seq_cov = col_double()
  )
)

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

fig_s1b_left <- swiss_model_stat %>% 
  ggplot(
    mapping = aes(
      x = "", 
      y = seq_id
    )
  ) +
  geom_violin(
    trim = FALSE,
    size = 1,
    fill = pnw_palette(name = "Shuksan2", n = 7)[c(4)]
  ) +
  geom_boxplot(
    width = 0.1,
    fill = "white",
    outlier.shape = NA
  ) +
  labs(
    x = "",
    y = str_wrap("Sequence identity (%)", width = 30)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 100), 
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s1", "fig_s1b_left.svg", sep = "/"),
  plot = fig_s1b_left,
  width = 4,
  height = 6,
  units = "in",
  device = "svg",
)


fig_s1b_right <- swiss_model_stat %>% 
  ggplot(
    mapping = aes(
      x = "", 
      y = qnorm
    )
  ) +
  geom_violin(
    trim = FALSE,
    size = 1,
    fill = pnw_palette(name = "Shuksan2", n = 7)[c(4)]
  ) +
  geom_boxplot(
    width = 0.1,
    fill = "white",
    outlier.shape = NA
  ) +
  labs(
    x = "",
    y = str_wrap("Normalized QMEAN", width = 30)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1),
    labels = sprintf("%.1f", seq(0, 1, 0.2)),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s1", "fig_s1b_right.svg", sep = "/"),
  plot = fig_s1b_right,
  width = 4,
  height = 6,
  units = "in",
  device = "svg",
)
