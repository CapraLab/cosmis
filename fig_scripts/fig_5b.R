# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())

# data and figure paths
data_dir <- "/path/to/datasets"
figure_dir <- "/path/to/where/figures/are/stored"

# load data sets
cosmis_df_pdb <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_pdb.tsv", sep = "/"), 
  col_names = TRUE,
  col_types = cols(
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
    syn_var_sites = col_integer(),
    total_syn_sites = col_double(),
    mis_var_sites = col_integer(),
    total_mis_sites = col_double(),
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
)

# load swiss model dataset
cosmis_df_swiss_model <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_swiss_model.tsv", sep = "/"),
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    seq_separations = col_character(),
    num_contacts = col_integer(),
    syn_var_sites = col_integer(),
    total_syn_sites = col_double(),
    mis_var_sites = col_integer(),
    total_mis_sites = col_double(),
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
)

# load alphafold dataset
cosmis_df_alphafold <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_alphafold.tsv", sep = "/"),
  col_names = TRUE,
  col_types = cols(
    uniprot_id = col_character(),
    enst_id = col_character(),
    uniprot_pos = col_integer(),
    uniprot_aa = col_character(),
    seq_separations = col_character(),
    num_contacts = col_integer(),
    syn_var_sites = col_integer(),
    total_syn_sites = col_double(),
    mis_var_sites = col_integer(),
    total_mis_sites = col_double(),
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
)

# combine the data set
comm_cols <- Reduce(
  intersect,
  list(
    colnames(cosmis_df_pdb), 
    colnames(cosmis_df_swiss_model),
    colnames(cosmis_df_alphafold)
  )
)
combined <- rbind(
  cosmis_df_pdb[comm_cols],
  cosmis_df_swiss_model[comm_cols],
  cosmis_df_alphafold[comm_cols]
)
combined <- combined %>% mutate(
  cosmis = (cs_mis_obs - mis_pmt_mean) / mis_pmt_sd,
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)

# load ClinVar variant IDs
clinvar_patho_ids <-  read_tsv(
  file = paste(
    result_dir, 
    "clinvar_unambiguous_pathogenic_ids_20210807.tsv", 
    sep = "/"
  ), 
  col_names = TRUE
)
clinvar_patho_ids <- clinvar_patho_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)
clinvar_benign_ids <- read_tsv(
  file = paste(
    result_dir, 
    "clinvar_unambiguous_benign_ids_20210807.tsv", 
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

# function for computing 
odds_ratio <- function(low, high) {
  if (low == 0 || high == 0) {
    a <- sum(patho$cosmis <= high & patho$cosmis >= low)
    c <- sum(benign$cosmis <= high & benign$cosmis >= low)
  } else {
    a <- sum(patho$cosmis <= high & patho$cosmis > low)
    c <- sum(benign$cosmis <= high & benign$cosmis > low)
  }
  b <- nrow(patho) - a
  d <- nrow(benign) - c
  or <- (a / b) / (c / d)
  ci_low <- exp(log(or) - 1.96 * sqrt(1 / a + 1 / b + 1 / c + 1 / d))
  ci_high <- exp(log(or) + 1.96 * sqrt(1 / a + 1 / b + 1 / c + 1 / d))
  return(c(or, ci_low, ci_high))
}

or_ci_m <- matrix(nrow = 10, ncol = 3)
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
  ci_high = or_ci_m[, 3]
)

# barplot theme
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
barplot_theme <- theme_classic() + theme(
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

# make a barplot
or_barplot <- ggplot(
  data = or_ci_df,
  mapping = aes(x = as.factor(percentile), y = or, fill = as.factor(percentile))
) + 
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(
    mapping = aes(ymin = ci_low, ymax = ci_high),
    width = 0.2
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "COSMIS score percentile",
    y = "OR (pathogenic vs benign)"
  ) +
  scale_y_continuous(
    limits = c(0, 15),
    breaks = seq(0, 15, 5),
    labels = sprintf("%.1f", seq(0, 15, 5)),
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
  barplot_theme

# save the bar plot to disk
ggsave(
  filename = paste(figure_dir, "clinvar_cosmis_or_barplot.svg", sep = "/"),
  plot = or_barplot,
  width = 8,
  height = 6,
  units = "in",
  device = "svg",
)