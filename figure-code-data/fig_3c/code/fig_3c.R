# load required packages
rm(list = ls())
library(tidyverse)
library(PNWColors)


# path to source data folder
data_dir <- "/path/to/source_data"

# load data set
fig_3c_data <- read_tsv(
  file = paste(data_dir, "fig_3c/data", "fig_3c_source_data.tsv", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    fraction = col_double()
  )
)

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

fig_3c <-  ggplot() + 
  geom_histogram(
    mapping = aes(fig_3c_data$fraction),
    binwidth = 0.05,
    boundary = 0.0,
    closed = "right",
    colour = "black",
    size = 0.25,
    fill = pnw_palette(name = "Shuksan2", n = 7)[2]
  ) +
  scale_x_continuous(
    name = str_wrap(
      "Fraction of high-confidence constrained sites", width = 30
    ),
    limits = c(0.0, 0.8),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Number of proteins",
    limits = c(0, 8000),
    breaks = seq(0, 8000, 2000),
    expand = c(0, 0)
  ) +
  hist_plot_theme

# save the histogram to disk
ggsave(
  filename = paste(data_dir, "fig_3c", "fig_3c.svg", sep = "/"),
  plot = fig_3c,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)
