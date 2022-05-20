rm(list = ls())
library(tidyverse)
library(PNWColors)


# path to source data folder
data_dir <- "/path/to/source_data"

enst_rate_count_df <- read_tsv(
  file = str_c(
    data_dir, "fig_s18/data", 
    "ensembl_cnl_enst_mp_and_counts.tsv",
    sep = "/"  
  ),
  col_names = TRUE
)

# fit a linear model
syn_lm_fit <- enst_rate_count_df %>% lm(
  formula = syn_count ~ syn_prob
)
summary(syn_lm_fit)

enst_rate_count_df <-  enst_rate_count_df %>% mutate(
  syn_exp = predict(
    object = syn_lm_fit, 
    newdata = data.frame(syn_prob = syn_prob)
  ),
  mis_exp = predict(
    object = syn_lm_fit, 
    newdata = data.frame(syn_prob = mis_prob)
  )
)

# make scatter plot
plot_margin <- margin(
  t = 0.5, r = 1.0, b = 0.5, l = 0.5, unit = "cm"
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

fig_s18a <- enst_rate_count_df %>% 
  filter(syn_exp < 4000) %>% ggplot(
    mapping = aes(
      x = syn_count,
      y = syn_exp
    )
  ) +
  geom_point(
    pch = 19,
    size = 1,
    color = pnw_palette(name = "Shuksan2", n = 5)[2]
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    size = 1,
    colour = "gray40"
  ) +
  scale_x_continuous(
    limits = c(0, 4000),
    breaks = seq(0, 4000, 1000),
    name = "",
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 4000),
    breaks = seq(0, 4000, 1000),
    name = "",
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(
    data_dir, "fig_s18", "fig_s18a.svg", sep = "/"
  ),
  plot = fig_s18a,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)

fig_s18b <- enst_rate_count_df %>% 
  filter(syn_exp < 4000) %>% ggplot(
    mapping = aes(
      x = mis_count,
      y = mis_exp
    )
  ) +
  geom_point(
    pch = 19,
    size = 1,
    color = pnw_palette(name = "Shuksan2", n = 5)[4]
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    size = 1,
    colour = "gray40"
  ) +
  scale_x_continuous(
    limits = c(0, 4000),
    breaks = seq(0, 4000, 1000),
    name = "",
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 4000),
    breaks = seq(0, 4000, 1000),
    name = "",
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(
    data_dir, "fig_s18", "fig_s18b.svg", sep = "/"
  ),
  plot = fig_s18b,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)
