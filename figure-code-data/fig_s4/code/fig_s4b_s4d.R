# load required packages
library(tidyverse)
library(PNWColors)
rm(list = ls())

################################################################################
# Depending on your computer, this script may take 5 hours to run!
################################################################################

# path to source data folder
data_dir <- "/path/to/source_data"

# load data sets
cosmis_df_pdb <- read_tsv(
  file = paste(data_dir, "cosmis_dataset_pdb_10a.tsv.gz", sep = "/"), 
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
  file = paste(data_dir, "cosmis_dataset_swiss_model_10a.tsv.gz", sep = "/"),
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
  file = paste(data_dir, "cosmis_dataset_alphafold_10a.tsv.gz", sep = "/"),
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

# function for computing fraction of long-range contacts
compute_lr_fraction <- function(x, cutoff = 15) {
  seq_dist <- str_split(x, pattern = ";") %>% unlist() %>% 
    parse_integer() %>% abs()
  # each contact was counted twice, once for A->B and once for B->A
  # here only retein one instance for each contact
  return(sum(seq_dist > cutoff) / length(seq_dist))
}

# compute fractios of long-range contacts
pdb_fracs <- cosmis_df_pdb %>% 
  .$seq_separations %>% 
  sapply(
    FUN = compute_lr_fraction,
    simplify = "array",
    USE.NAMES = FALSE
  )
swiss_model_fracs <- cosmis_df_swiss_model %>% 
  .$seq_separations %>% 
  sapply(
    FUN = compute_lr_fraction,
    simplify = "array",
    USE.NAMES = FALSE
  )
af2_fracs <- cosmis_df_alphafold %>% 
  .$seq_separations %>% 
  sapply(
    FUN = compute_lr_fraction,
    simplify = "array",
    USE.NAMES = FALSE
  )

long_range_contact_frac <- bind_rows(
  tibble(
    dataset = "PDB",
    fxn = pdb_fracs
  ),
  tibble(
    dataset = "SWISS-MODEL",
    fxn = swiss_model_fracs
  ),
  tibble(
    dataset = "AF2",
    fxn = af2_fracs
  )
)

lr_fxns <- data.frame(
  class = c(">10%", ">20%", ">30%", ">40%", ">50%"),
  fxn = c(
    mean(long_range_contact_frac$fxn >= 0.1, na.rm = TRUE),
    mean(long_range_contact_frac$fxn >= 0.2, na.rm = TRUE),
    mean(long_range_contact_frac$fxn >= 0.3, na.rm = TRUE),
    mean(long_range_contact_frac$fxn >= 0.4, na.rm = TRUE),
    mean(long_range_contact_frac$fxn >= 0.5, na.rm = TRUE)
  )
)

# make a bar plot
plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"
)
plot_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text = element_text(
    size = 20, color = "black"
  ),
  axis.title.x = element_text(
    color = "black", size = 24, margin = margin(t = 10)
  ),
  axis.title.y = element_text(
    color = "black", size = 24, margin = margin(r = 10)
  ),
  axis.ticks.length = unit(.25, "cm"),
  axis.ticks = element_line(),
  legend.position = "none",
  plot.margin = plot_margin
)

fig_s4b <- ggplot(
  data = lr_fxns, 
  mapping = aes(
    x = class, 
    y = fxn,
    fill = class
  )
) +
  geom_bar(
    stat = "identity"
  ) +
  scale_fill_manual(
    values = rep(pnw_palette(name = "Shuksan2", n = 5)[2], 5)
  ) + 
  xlab(
    label = str_wrap(
      "Fraction of contacts a site makes that are long-range",
      width = 35
    )
  ) +
  scale_y_continuous(
    name = str_wrap("Fraction of all 6.1 million sites", width = 40),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s4", "fig_s4b.svg", sep = "/"),
  plot = fig_s4b,
  width = 6,
  height = 6,
  units = "in",
  device = "svg",
)


################################################################################
# Compute the fraction of neighboring 30 residues (15 on each side) 
# along the sequence that are not in contact with the central site.
################################################################################

# function for computing fractions of sequence neighbors 
# that are not in 3D contact
fxn_not_in_contact <- function(x, cutoff = 15) {
  seq_dist <- str_split(x, pattern = ";") %>% unlist %>% 
    parse_integer() %>% abs()
  return(1 - sum(seq_dist <= cutoff) / (cutoff * 2))
}

pdb_fxn_not_contact <- sapply(
  X = cosmis_df_pdb$seq_separations, 
  FUN = fxn_not_in_contact, 
  simplify = "array", 
  USE.NAMES = FALSE
)
sm_fxn_not_contact <- sapply(
  X = cosmis_df_swiss_model$seq_separations, 
  FUN = fxn_not_in_contact, 
  simplify = "array", 
  USE.NAMES = FALSE
)
af2_fxn_not_in_contact <- sapply(
  X = cosmis_df_alphafold$seq_separations, 
  FUN = fxn_not_in_contact, 
  simplify = "array", 
  USE.NAMES = FALSE
)

fxns_nic <- bind_rows(
  tibble(
    dataset = "PDB",
    fxn = pdb_fxn_not_contact
  ),
  tibble(
    dataset = "SWISS-MODEL",
    fxn = sm_fxn_not_contact
  ),
  tibble(
    dataset = "AF2",
    fxn = af2_fxn_not_in_contact
  )
)

nic_fxns_df <- data.frame(
  class = c(">50%", ">60%", ">70%", ">80%", ">90%"),
  fxn = c(
    mean(fxns_nic$fxn >= 0.5, na.rm = TRUE),
    mean(fxns_nic$fxn >= 0.6, na.rm = TRUE),
    mean(fxns_nic$fxn >= 0.7, na.rm = TRUE),
    mean(fxns_nic$fxn >= 0.8, na.rm = TRUE),
    mean(fxns_nic$fxn >= 0.9, na.rm = TRUE)
  )
)

fig_s4d <- ggplot(
  data = nic_fxns_df, 
  mapping = aes(
    x = class, 
    y = fxn,
    fill = class
  )
) +
  geom_bar(
    stat = "identity"
  ) +
  scale_fill_manual(
    values = rep(pnw_palette(name = "Shuksan2", n = 5)[2], 5)
  ) + 
  xlab(
    label = str_wrap(
      "Fraction of neighboring sites that are not in 3D contact",
      width = 35
    )
  ) +
  scale_y_continuous(
    name = str_wrap("Fraction of all 6.1 million sites", width = 40),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  plot_theme

ggsave(
  filename = paste(data_dir, "fig_s4", "fig_s4d.svg", sep = "/"),
  plot = fig_s4d,
  width = 6,
  height = 6,
  units = "in",
  device = "svg",
)
