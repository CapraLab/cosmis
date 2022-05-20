# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())

# path to source data folder
data_dir <- "/path/to/source_data"

# load data sets
cosmis_pdb <- read_tsv(
  file = paste(data_dir, "cosmis_scores_pdb.tsv.gz", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    cosmis = col_double(),
    cosmis_pvalue = col_double(),
    cossyn = col_double(),
    str_source = col_character()
  )
)
cosmis_swiss_model <- read_tsv(
  file = paste(data_dir, "cosmis_scores_swiss_model.tsv.gz", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    cosmis = col_double(),
    cosmis_pvalue = col_double(),
    cossyn = col_double(),
    str_source = col_character()
  )
)
cosmis_alphafold <- read_tsv(
  file = paste(data_dir, "cosmis_scores_alphafold.tsv.gz", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    cosmis = col_double(),
    cosmis_pvalue = col_double(),
    cossyn = col_double(),
    str_source = col_character()
  )
)
combined <- bind_rows(
  cosmis_pdb,
  cosmis_swiss_model,
  cosmis_alphafold
)
combined <- combined %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)

# load ClinVar variant IDs
clinvar_patho_ids <-  read_tsv(
  file = paste(
    data_dir, "fig_5b/data", 
    "clinvar_unambiguous_all_pathogenic_ids_20210807.tsv", 
    sep = "/"
  ), 
  col_names = TRUE
)
clinvar_patho_ids <- clinvar_patho_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)
clinvar_benign_ids <- read_tsv(
  file = paste(
    data_dir, "fig_5b/data",
    "clinvar_unambiguous_all_benign_ids_20210807.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)
clinvar_benign_ids <- clinvar_benign_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)

# cosmis scores for ClinVar pathogenic variants
patho_cosmis <- combined %>% filter(
  full_id %in% clinvar_patho_ids$full_id
) %>% mutate(
  label = 1
)
benign_cosmis <- combined %>% filter(
  full_id %in% clinvar_benign_ids$full_id
) %>% mutate(
  label = 0
)

cosmis_clinvar <- rbind(
  patho_cosmis,
  benign_cosmis
)
cosmis_clinvar <- cosmis_clinvar %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )


# COSMIS score percentiles
percentiles <- quantile(combined$cosmis, seq(0.1, 1, 0.1))

# function for computing odds ratio
odds_ratio <- function(low, high) {
  if (low == 0 || high == 0) {
    a <- sum(patho_cosmis$cosmis <= high & patho_cosmis$cosmis >= low)
    c <- sum(benign_cosmis$cosmis <= high & benign_cosmis$cosmis >= low)
  } else {
    a <- sum(patho_cosmis$cosmis <= high & patho_cosmis$cosmis > low)
    c <- sum(benign_cosmis$cosmis <= high & benign_cosmis$cosmis > low)
  }
  b <- nrow(patho_cosmis) - a
  d <- nrow(benign_cosmis) - c
  or <- (a / b) / (c / d)
  ci_low <- exp(log(or) - 1.96 * sqrt(1 / a + 1 / b + 1 / c + 1 / d))
  ci_high <- exp(log(or) + 1.96 * sqrt(1 / a + 1 / b + 1 / c + 1 / d))
  return(c(or, ci_low, ci_high, a, b, c, d))
}

or_ci_m <- matrix(nrow = 10, ncol = 7)
for (i in 1:9) {
  or_ci <- odds_ratio(low = percentiles[i], high = percentiles[i + 1])
  or_ci_m[i + 1, ] <- or_ci
}
# most constrained percentile
x <- odds_ratio(low = -6, high = percentiles[1])
or_ci_m[1, ] <- x

# create a dataframe for plotting
or_ci_df <- data.frame(
  percentile = seq(0.1, 1, 0.1),
  or = or_ci_m[, 1],
  ci_low = or_ci_m[, 2],
  ci_high = or_ci_m[, 3],
  a = or_ci_m[, 4],
  b = or_ci_m[, 5],
  c = or_ci_m[, 6],
  d = or_ci_m[, 7]
)
write_csv(
  or_ci_df, 
  file = str_c(data_dir, "fig_5b", "odds_ratios.csv", sep = "/")
)

# barplot theme
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme_classic() + theme(
  axis.text = element_text(
    size = 16, color = "black"
  ),
  axis.text.x.bottom = element_text(
    angle = 90, vjust = 0.5, hjust = 1
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

fig_5b <- ggplot(
  data = or_ci_df,
  mapping = aes(x = as.factor(percentile), y = or, fill = as.factor(percentile))
) +
  geom_errorbar(
    mapping = aes(ymin = ci_low, ymax = ci_high),
    width = 0.2
  ) +
  geom_point(size = 5, shape = 21) +

  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "COSMIS score percentile",
    y = "OR (pathogenic vs benign)"
  ) +
  scale_y_continuous(
    limits = c(0, 12),
    breaks = seq(0, 12, 3),
    labels = sprintf("%.1f", seq(0, 12, 3)),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    labels = c(
      " 0 - 10",
      "10 - 20",
      "20 - 30",
      "30 - 40",
      "40 - 50",
      "50 - 60",
      "60 - 70",
      "70 - 80",
      "80 - 90",
      "90 - 100"
    )
  ) +
  scale_fill_manual(
    values = rep(pnw_palette(name = "Shuksan2", n = 7)[2], 10)
  ) + 
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_5b", "fig_5b.svg", sep = "/"),
  plot = fig_5b,
  width = 10,
  height = 6,
  units = "in",
  device = "svg",
)