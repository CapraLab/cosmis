# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())


# data and figure paths
data_dir <- "/path/to/datasets/"
figure_dir <- "/path/to/where/figures/are/stored/"

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
  full_id = str_c(uniprot_id, uniprot_pos, sep = "")
)

#===============================================================================
# load other scores
#===============================================================================
# gerp score
pdb_gerp <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_gerp_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_gerp <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_gerp_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_gerp <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_gerp_scores.tsv", sep = "/"),
  col_names = TRUE
)
gerp_scores <- bind_rows(
  pdb_gerp, swiss_model_gerp, alphafold_gerp
)

# phyloP scores
pdb_phylop <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_phylop_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_phylop <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_phylop_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_phylop <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_phylop_scores.tsv", sep = "/"),
  col_names = TRUE
)
phylop_scores <- bind_rows(
  pdb_phylop, swiss_model_phylop, alphafold_phylop
)

# phastCons scores
pdb_phastconst <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_phastcons_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_phastconst <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_phastcons_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_phastconst <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_phastcons_scores.tsv", sep = "/"),
  col_names = TRUE
)
phastcons_scores <- bind_rows(
  pdb_phastconst, swiss_model_phastconst, alphafold_phastconst
)

# ConSurf scores
pdb_consurf <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_consurf_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_consurf <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_consurf_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_consurf <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_consurf_scores.tsv", sep = "/"),
  col_names = TRUE
)
consurf_scores <- bind_rows(
  pdb_consurf, swiss_model_consurf, alphafold_consurf
)

# pLI missense_z
pdb_pli <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_pli_mis_z_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_pli <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_pli_mis_z_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_pli <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_pli_mis_z_scores.tsv", sep = "/"),
  col_names = TRUE
)
pli_scores <- bind_rows(
  pdb_pli, swiss_model_pli, alphafold_pli
)

# RVIS scores
pdb_rvis <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_rvis_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_rvis <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_rvis_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_rvis <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_rvis_scores.tsv", sep = "/"),
  col_names = TRUE,
  col_types = cols(
    hgnc = col_character(),
    rvis = col_double(),
    rvis_percentile = col_double()
  )
)
rvis_scores <- bind_rows(
  pdb_rvis, swiss_model_rvis, alphafold_rvis
)

# MTR scores
pdb_mtr <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_mtr_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_mtr <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_mtr_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_mtr <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_mtr_scores.tsv", sep = "/"),
  col_names = TRUE
)
mtr_scores <- bind_rows(
  pdb_mtr, swiss_model_mtr, alphafold_mtr
)

# MTR3D scores
pdb_mtr3d <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_pdb_mtr3d_scores.tsv", sep = "/"),
  col_names = TRUE
)
swiss_model_mtr3d <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_swiss_model_mtr3d_scores.tsv", sep = "/"),
  col_names = TRUE
)
alphafold_mtr3d <- read_tsv(
  file = str_c(data_dir, "cosmis_dataset_alphafold_mtr3d_scores.tsv", sep = "/"),
  col_names = TRUE
)
mtr3d_scores <- bind_rows(
  pdb_mtr3d, swiss_model_mtr3d, alphafold_mtr3d
)


# combine the data set
select <- dplyr::select
combined_na_removed <- bind_cols(
  combined %>% select(uniprot_id, enst_id, uniprot_pos, uniprot_aa, cosmis), gerp_scores %>% select(gerp),
  phylop_scores %>% select(phylop), phastcons_scores %>% select(phastcons),
  consurf_scores %>% select(consurf), pli_scores %>% select(pli, mis_z),
  mtr_scores %>% select(mtr), rvis_scores %>% select(rvis), mtr3d_scores %>% select(mtr3d)
) %>% drop_na()

# load de novo variants coordinates
case_coords <- read_table(
  file = str_c(
    data_dir, "de_novo_dataset_samocha_case_variant_coord.txt", sep = "/"
  ),
  col_names = TRUE
)
control_coords <- read_table(
  file = str_c(
    data_dir, "de_novo_dataset_samocha_control_variant_coord.txt", sep = "/"
  ),
  col_names = TRUE
)

combined_na_removed <- combined_na_removed %>% mutate(
  keys = str_c(enst_id, uniprot_pos, sep = "_")
)

control_cosmis_df <- combined_na_removed %>% filter(
  keys %in% control_coords$coord
)
case_cosmis_df <- combined_na_removed %>% filter(
  keys %in% case_coords$coord
)

case_cosmis_df <- case_cosmis_df %>% mutate(
  class = "case"
)
control_cosmis_df <- control_cosmis_df %>% mutate(
  class = "control"
)

de_novo_cosmis <- bind_rows(
  case_cosmis_df, control_cosmis_df
) %>% dplyr::select(-keys)


# define a function for computing odds ratios
odds_ratio <- function(score, threshold, direction = "low") {
  if(direction == "low") {
    a <- sum(case_cosmis_df %>% select(score) <= threshold)
    b <- sum(control_cosmis_df %>% select(score) <= threshold)
    c <- sum(case_cosmis_df %>% select(score) > threshold)
    d <- sum(control_cosmis_df %>% select(score) > threshold)
  } else {
    a <- sum(case_cosmis_df %>% select(score) >= threshold)
    b <- sum(control_cosmis_df %>% select(score) >= threshold)
    c <- sum(case_cosmis_df %>% select(score) < threshold)
    d <- sum(control_cosmis_df %>% select(score) < threshold)    
  }
  dat <- data.frame(
    pathogenic = c(a, b),
    benign = c(c, d),
    row.names = c("Case", "Control"),
    stringsAsFactors = FALSE
  )
  f_test <- fisher.test(dat)
  print(f_test)
  return(c(f_test$estimate, f_test$p.value, f_test$conf.int))
}

# compute percentiles
cosmis_percentiles <- quantile(combined$cosmis, seq(0.05, 1, 0.05))
mtr_percentiles <- quantile(mtr_scores$mtr, seq(0.05, 1, 0.05), na.rm = TRUE)
mtr3d_percentiles <- quantile(mtr3d_scores$mtr3d, seq(0.05, 1, 0.05), na.rm = TRUE)
rvis_percentiles <- quantile(rvis_scores$rvis, seq(0.05, 1, 0.05), na.rm = TRUE)
consurf_percentiles <- quantile(consurf_scores$consurf, seq(0.05, 1, 0.05), na.rm = TRUE)
pli_percentiles <- quantile(pli_scores$pli, seq(0.05, 1, 0.05), na.rm = TRUE)
mis_z_percentiles <- quantile(pli_scores$mis_z, seq(0.05, 1, 0.05), na.rm = TRUE)
gerp_percentiles <- quantile(gerp_scores$gerp, seq(0.05, 1, 0.05), na.rm = TRUE)
phylop_percentiles <- quantile(phylop_scores$phylop, seq(0.05, 1, 0.05), na.rm = TRUE)
phastcons_percentiles <- quantile(phastcons_scores$phastcons, seq(0.05, 1, 0.05), na.rm = TRUE)

# now compute odds ratios
cosmis_or <- odds_ratio("cosmis", cosmis_percentiles[1])
consurf_or <- odds_ratio("consurf", consurf_percentiles[1])
mtr_or <- odds_ratio("mtr", mtr_percentiles[1])
mtr3d_or <- odds_ratio("mtr3d", mtr3d_percentiles[1])
rvis_or <- odds_ratio("rvis", rvis_percentiles[1])
pli_or <- odds_ratio("pli", tail(pli_percentiles, 5)[4], direction = "high")
mis_z_or <- odds_ratio("mis_z", tail(mis_z_percentiles, 5)[4], direction = "high")
gerp_or <- odds_ratio("gerp", tail(gerp_percentiles, 5)[4], direction = "high")
phylop_or <- odds_ratio("phylop", tail(phylop_percentiles, 5)[4], direction = "high")
phastcons_or <- odds_ratio("phastcons", tail(phastcons_percentiles, 5)[4], direction = "high")

# combine odds ratios
all_or <- rbind(
  cosmis_or,
  consurf_or,
  mtr_or,
  mtr3d_or,
  rvis_or,
  pli_or,
  mis_z_or,
  gerp_or,
  phylop_or,
  phastcons_or
)
colnames(all_or) <-   c("odds_ratio", "p_value", "low", "high")
all_or <- as.tibble(all_or) %>% mutate(
  metric = c(
    "COSMIS", "ConSurf", "MTR", "MTR3D", "RVIS", "pLI", "Missense Z", 
    "GERP++", "phyloP", "phastCons"
  )
)

# barplot theme
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
barplot_theme <- theme_classic() + theme(
  axis.text.x = element_text(
    size = 16, color = "black", angle = 90, hjust = 1, vjust = 0.5
  ),
  axis.text.y = element_text(
    size = 16, color = "black"
  ),
  axis.title.x = element_blank(),
  axis.title.y = element_text(
    color = "black", size = 20, margin = margin(r = 10)
  ),
  axis.ticks.length = unit(.25, "cm"),
  axis.ticks = element_line(),
  legend.position = "none",
  plot.margin = plot_margin
)

# bin the cosmis score
or_barplot <- all_or %>% arrange(desc(odds_ratio)) %>% 
  mutate(
    metric = factor(metric, levels = metric)
  ) %>%
  ggplot(
    mapping = aes(
      x = metric, 
      y = odds_ratio
    )
  ) + 
  geom_bar(
    stat = "identity", 
    width = 0.7,
    fill = pnw_palette(name = "Shuksan2", n = 7)[2]
  ) +
  geom_errorbar(
    mapping = aes(ymin = low, ymax = high),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  # geom_hline(yintercept = 1) +
  labs(
    y = "Odds ratio"
  ) +
  scale_y_continuous(
    limits = c(0, 8),
    breaks = seq(0, 8, 2),
    labels = sprintf("%.1f", seq(0, 8, 2)),
    expand = c(0, 0)
  ) +
  barplot_theme

# save the scatter plot to disk
ggsave(
  filename = paste(figure_dir, "de_novo_cosmis_or_barplot_5_pctl.svg", sep = "/"),
  plot = or_barplot,
  width = 10,
  height = 6,
  units = "in",
  device = "svg",
)
