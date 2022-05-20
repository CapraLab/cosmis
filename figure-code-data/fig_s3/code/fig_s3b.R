# load required packages
rm(list = ls())
library(tidyverse)
library(PNWColors)


# path to source data folder
data_dir <- "/path/to/source_data"

# load data sets
cosmis_df_pdb <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_pdb.tsv.gz", sep = "/"), 
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
  file = paste(data_dir, "cosmis_dataset_swiss_model.tsv.gz", sep = "/"),
  col_names = TRUE,
  col_types = cols(
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
)

# load alphafold dataset
cosmis_df_alphafold <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_alphafold.tsv.gz", sep = "/"),
  col_names = TRUE,
  col_types = cols(
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
)

count_global_contacts <- function(x, cutoff = 15) {
  seq_dist <- str_split(x, pattern = ";") %>% unlist() %>% 
    parse_integer() %>% abs()
  return(sum(seq_dist > cutoff))
}

pdb_gc_site_fractions <- cosmis_df_pdb %>% 
  mutate(
    num_gc = sapply(
      seq_separations,
      FUN = count_global_contacts,
      simplify = "array",
      USE.NAMES = FALSE
    )
  ) %>% 
  group_by(uniprot_id) %>% 
  summarise(frac = sum(num_gc >= 1, na.rm = TRUE) / n()) %>% 
  filter(frac > 0)

sm_gc_site_fraction <- cosmis_df_swiss_model %>% 
  mutate(
    num_gc = sapply(
      seq_separations,
      FUN = count_global_contacts,
      simplify = TRUE,
      USE.NAMES = FALSE
    )
  ) %>% 
  group_by(uniprot_id) %>% 
  summarise(frac = mean(num_gc >= 1, na.rm = TRUE)) %>% 
  filter(frac > 0)

af2_gc_site_fraction <- cosmis_df_alphafold %>% 
  mutate(
    num_gc = sapply(
      seq_separations,
      FUN = count_global_contacts,
      simplify = TRUE,
      USE.NAMES = FALSE
    )
  ) %>% 
  group_by(uniprot_id) %>% 
  summarise(frac = mean(num_gc >= 1, na.rm = TRUE)) %>% 
  filter(frac > 0)


per_prot_gc_fraction <- bind_rows(
  tibble(
    fraction = pdb_gc_site_fractions$frac,
    dataset = "PDB"
  ),
  tibble(
    fraction = sm_gc_site_fraction$frac,
    dataset = "SWISS-MODEL"
  ),
  tibble(
    fraction = af2_gc_site_fraction$frac,
    dataset = "AF2"
  )
)

plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme_classic() + theme(
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

fig_3b_pdb <- ggplot() + 
  geom_histogram(
    mapping = aes(pdb_gc_site_fractions$frac),
    binwidth = 0.025,
    boundary = 0,
    colour = "black",
    fill = pnw_palette(name = "Shuksan2", n = 7)[2]
  ) +
  scale_x_continuous(
    name = str_wrap(
      "Fraction of residues with long-range contact", width = 30
    ),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Number of proteins",
    limits = c(0, 1000),
    breaks = seq(0, 1000, 200),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s3", "fig_3b_pdb.svg", sep = "/"),
  plot = fig_3b_pdb,
  width = 6,
  height = 6,
  units = "in",
  device = "svg",
)

fig_3b_swiss_model <- ggplot() + 
  geom_histogram(
    mapping = aes(sm_gc_site_fraction$frac),
    binwidth = 0.025,
    boundary = 0,
    colour = "black",
    fill = pnw_palette(name = "Shuksan2", n = 7)[4]
  ) +
  scale_x_continuous(
    name = str_wrap(
      "Fraction of residues with long-range contact", width = 30
    ),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Number of proteins",
    limits = c(0, 1000),
    breaks = seq(0, 1000, 200),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s3", "fig_3b_swiss_model.svg", sep = "/"),
  plot = fig_3b_swiss_model,
  width = 6,
  height = 6,
  units = "in",
  device = "svg",
)

fig_3b_alphafold <- ggplot() + 
  geom_histogram(
    mapping = aes(af2_gc_site_fraction$frac),
    binwidth = 0.025,
    boundary = 0,
    colour = "black",
    fill = pnw_palette(name = "Shuksan2", n = 7)[6]
  ) +
  scale_x_continuous(
    name = str_wrap(
      "Fraction of residues with long-range contact", width = 30
    ),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Number of proteins",
    limits = c(0, 1000),
    breaks = seq(0, 1000, 200),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s3", "fig_3b_alphafold.svg", sep = "/"),
  plot = fig_3b_alphafold,
  width = 6,
  height = 6,
  units = "in",
  device = "svg",
)