# load required libraries
rm(list = ls())
library(plotROC)
library(tidyverse)
library(PNWColors)


# path to source data folder
data_dir <- "/path/to/source_data/"

# load ConSurf scores
patho_consurf <- read_tsv(
  file = paste(
    data_dir, "fig_s12/data",
    "clinvar_unambiguous_pathogenic_consurf_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_consurf$label <- 1
benign_consurf <- read_tsv(
  file = paste(
    data_dir, "fig_s12/data",
    "clinvar_unambiguous_benign_consurf_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_consurf$label <- 0
consurf_clinvar <- bind_rows(
  patho_consurf,
  benign_consurf
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load RSA scores
patho_rsa <- read_tsv(
  file = paste(
    data_dir, "fig_s12/data",
    "clinvar_unambiguous_pathogenic_rsa_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE,
  col_types = cols(
    rsa = col_double()
  )
)
patho_rsa$label <- 1
benign_rsa <- read_tsv(
  file = paste(
    data_dir, "fig_s12/data",
    "clinvar_unambiguous_benign_rsa_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE,
  col_types = cols(
    rsa = col_double()
  )
)
benign_rsa$label <- 0
rsa_clinvar <- bind_rows(
  patho_rsa,
  benign_rsa
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load phyloP scores
patho_phylop <- read_tsv(
  file = paste(
    data_dir, "fig_s12/data",
    "clinvar_unambiguous_pathogenic_phylop_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_phylop$label <- 1
benign_phylop <- read_tsv(
  file = paste(
    data_dir, "fig_s12/data",
    "clinvar_unambiguous_benign_phylop_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_phylop$label <- 0
phylop_clinvar <- bind_rows(
  patho_phylop,
  benign_phylop
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load phastCons scores
patho_phastcons <- read_tsv(
  file = paste(
    data_dir, "fig_s12/data",
    "clinvar_unambiguous_pathogenic_phastcons_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_phastcons$label <- 1
benign_phastcons <- read_tsv(
  file = paste(
    data_dir, "fig_s12/data",
    "clinvar_unambiguous_benign_phastcons_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_phastcons$label <- 0
phastcons_clinvar <- bind_rows(
  patho_phastcons,
  benign_phastcons
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load GERP scores
patho_gerp <- read_tsv(
  file = paste(
    data_dir, "fig_s12/data",
    "clinvar_unambiguous_pathogenic_gerp_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_gerp$label <- 1
benign_gerp <- read_tsv(
  file = paste(
    data_dir, "fig_s12/data",
    "clinvar_unambiguous_benign_gerp_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_gerp$label <- 0
gerp_clinvar <- bind_rows(
  patho_gerp,
  benign_gerp
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

cons_scores <- bind_cols(
  consurf_clinvar %>% select(class), consurf_clinvar %>% select(consurf),
  rsa_clinvar %>% select(rsa), phylop_clinvar %>% select(phylop),
  gerp_clinvar %>% select(gerp), phastcons_clinvar %>% select(phastcons)
) %>% drop_na() 
positives <- cons_scores %>% filter(class == "pathogenic")
negatives <- cons_scores %>% filter(class == "benign")

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
  axis.ticks = element_line(size = 1.0),
  legend.position = "none",
  aspect.ratio = 1.0,
  plot.margin = plot_margin
)

# ROC curves
fig_s12 <- ggplot() +
  geom_roc(
    data = cons_scores,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = - consurf
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = "blue"
  ) +
  geom_roc(
    data = cons_scores,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = gerp
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 11)[3]
  ) +
  geom_roc(
    data = cons_scores,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = phylop
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 11)[5]
  ) +
  geom_roc(
    data = cons_scores,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = phastcons
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 11)[7]
  ) +
  geom_roc(
    data = cons_scores,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = -rsa
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 11)[9]
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

# save the scatter plot to disk
ggsave(
  filename = paste(data_dir, "fig_s12", "fig_s12.svg", sep = "/"),
  plot = fig_s12,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)