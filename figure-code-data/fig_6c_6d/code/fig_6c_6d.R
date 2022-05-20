# load required libraries
rm(list = ls())
library(plotROC)
library(tidyverse)
library(PNWColors)
library(ROCR)


# path to source data folder
data_dir <- "/path/to/source_data/"

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
    data_dir, "fig_6c_6d/data",
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
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_all_benign_ids_20210807.tsv", 
    sep = "/"
  ),
  col_names = TRUE
)
clinvar_benign_ids <- clinvar_benign_ids %>% mutate(
  full_id = str_c(uniprot_id, uniprot_pos, sep = "_")
)
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

# combine ClinVar datasets
cosmis_clinvar <- bind_rows(
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

# load pLI and Missense_Z scores
patho_pli_mis_z <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_pathogenic_pli_mis_z_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_pli_mis_z$label <- 1
benign_pli_mis_z <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_benign_pli_mis_z_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_pli_mis_z$label <- 0
pli_mis_z_clinvar <- bind_rows(
  patho_pli_mis_z,
  benign_pli_mis_z
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load MTR3D scores
patho_mtr3d <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_pathogenic_mtr3d_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_mtr3d$label <- 1
benign_mtr3d <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_benign_mtr3d_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_mtr3d$label <- 0
mtr3d_clinvar <- bind_rows(
  patho_mtr3d,
  benign_mtr3d
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load MTR scores
patho_mtr <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_pathogenic_mtr_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_mtr$label <- 1
benign_mtr <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_benign_mtr_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_mtr$label <- 0
mtr_clinvar <- bind_rows(
  patho_mtr,
  benign_mtr
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load RVIS scores
patho_rvis <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_pathogenic_rvis_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_rvis$label <- 1
benign_rvis <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_benign_rvis_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
benign_rvis$label <- 0
rvis_clinvar <- bind_rows(
  patho_rvis,
  benign_rvis
) %>% 
  mutate(
    class = factor(
      ifelse(label == 0, "benign", "pathogenic"), 
      levels = c("benign", "pathogenic")
    )
  )

# load ConSurf scores
patho_consurf <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_pathogenic_consurf_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_consurf$label <- 1
benign_consurf <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
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
    data_dir, "fig_6c_6d/data",
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
    data_dir, "fig_6c_6d/data",
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
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_pathogenic_phylop_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_phylop$label <- 1
benign_phylop <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
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
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_pathogenic_phastcons_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_phastcons$label <- 1
benign_phastcons <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
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
    data_dir, "fig_6c_6d/data",
    "clinvar_unambiguous_pathogenic_gerp_20210807.tsv",
    sep = "/"
  ),
  col_names = TRUE
)
patho_gerp$label <- 1
benign_gerp <- read_tsv(
  file = paste(
    data_dir, "fig_6c_6d/data",
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

# COSMIS vs other metrics
cosmis_vs_others <- bind_cols(
  cosmis_clinvar, mtr3d_clinvar %>% select(mtr3d),
  rvis_clinvar %>% select(rvis), pli_mis_z_clinvar %>% select(pli, mis_z),
  mtr_clinvar %>% select(mtr), consurf_clinvar %>% select(consurf),
  rsa_clinvar %>% select(rsa), phylop_clinvar %>% select(phylop),
  gerp_clinvar %>% select(gerp), phastcons_clinvar %>% select(phastcons)
) %>% drop_na()
positives <- cosmis_vs_others %>% filter(class == "pathogenic")
negatives <- cosmis_vs_others %>% filter(class == "benign")

# repeated k-fold cross-validation
r <- 1 # number of repeats
k <- 5 # fold of cross-validation
aucs_model1 <- numeric(length = k)
cv_results_model1 <- data.frame(cv_fold = character(), pred = double(), class = character())
set.seed(seed = 42)
for(i in 1:r) {
  # shuffle the subsets
  negatives <- negatives[sample(1:nrow(negatives)),]
  positives <- positives[sample(1:nrow(positives)),]
  neg.folds <- cut(x = 1:nrow(negatives), breaks = k, labels = FALSE)
  pos.folds <- cut(x = 1:nrow(positives), breaks = k, labels = FALSE)
  for(j in 1:k) {
    # test set
    test <- rbind(negatives[neg.folds == j,], positives[pos.folds == j,])
    # training set
    train <- rbind(negatives[neg.folds !=j,], positives[pos.folds != j,])
    # logistic regression model
    logit.model <- glm(
      class ~ rvis + pli + mis_z + mtr + mtr3d, 
      data = train, family = "binomial"
    )
    # make predictions on the test set
    test$logit.pred <- predict(object = logit.model, newdata = test, type = "response")
    # compute area under the ROC curve, or AUC
    rocr <- prediction(predictions = test$logit.pred, labels = test$class)
    perf <- performance(prediction.obj = rocr, measure = "auc")
    aucs_model1[(i-1)*k+j] <- unlist(perf@y.values)
    # save cv results
    tmp_df <- data.frame(
      cv_fold = paste("fold", j, sep = ""), 
      pred = test$logit.pred,
      class = test$class
    )
    cv_results_model1 <- rbind(cv_results_model1, tmp_df)
  }
}

# repeated k-fold cross-validation
r <- 1
k <- 5
aucs_model2 <- numeric(length = k)
cv_results_model2 <- data.frame(
  cv_fold = character(), pred = double(), class = character()
)
set.seed(seed = 42)
for(i in 1:r) {
  # shuffle the subsets
  negatives <- negatives[sample(1:nrow(negatives)),]
  positives <- positives[sample(1:nrow(positives)),]
  neg.folds <- cut(x = 1:nrow(negatives), breaks = k, labels = FALSE)
  pos.folds <- cut(x = 1:nrow(positives), breaks = k, labels = FALSE)
  for(j in 1:k) {
    # test set
    test <- rbind(negatives[neg.folds == j,], positives[pos.folds == j,])
    # training set
    train <- rbind(negatives[neg.folds !=j,], positives[pos.folds != j,])
    # logistic regression model
    logit.model <- glm(
      # class ~ rvis + pli + mis_z + mtr + mtr3d,
      class ~ consurf + cosmis,
      data = train, family = "binomial"
    )
    # make predictions on the test set
    test$logit.pred <- predict(object = logit.model, newdata = test, type = "response")
    # compute area under the ROC curve, or AUC
    rocr <- prediction(predictions = test$logit.pred, labels = test$class)
    perf <- performance(prediction.obj = rocr, measure = "auc")
    aucs_model2[(i-1)*k+j] <- unlist(perf@y.values)
    # save cv results
    tmp_df <- data.frame(
      cv_fold = paste("fold", j, sep = ""), 
      pred = test$logit.pred,
      class = test$class
    )
    cv_results_model2 <- rbind(cv_results_model2, tmp_df)
  }
}

# repeated k-fold cross-validation
r <- 1
k <- 5
aucs_model3 <- numeric(length = k)
cv_results_model3 <- data.frame(cv_fold = character(), pred = double(), class = character())
set.seed(seed = 42)
for(i in 1:r) {
  # shuffle the subsets
  negatives <- negatives[sample(1:nrow(negatives)),]
  positives <- positives[sample(1:nrow(positives)),]
  neg.folds <- cut(x = 1:nrow(negatives), breaks = k, labels = FALSE)
  pos.folds <- cut(x = 1:nrow(positives), breaks = k, labels = FALSE)
  for(j in 1:k) {
    # test set
    test <- rbind(negatives[neg.folds == j,], positives[pos.folds == j,])
    # training set
    train <- rbind(negatives[neg.folds !=j,], positives[pos.folds != j,])
    # logistic regression model
    logit.model <- glm(
      class ~ cosmis + rvis + pli + mis_z + mtr + mtr3d + consurf + phylop + gerp + phastcons + rsa, 
      data = train, family = "binomial"
    )
    # make predictions on the test set
    test$logit.pred <- predict(object = logit.model, newdata = test, type = "response")
    # compute area under the ROC curve, or AUC
    rocr <- prediction(predictions = test$logit.pred, labels = test$class)
    perf <- performance(prediction.obj = rocr, measure = "auc")
    aucs_model3[(i-1)*k+j] <- unlist(perf@y.values)
    # save cv results
    tmp_df <- data.frame(
      cv_fold = paste("fold", j, sep = ""), 
      pred = test$logit.pred,
      class = test$class
    )
    cv_results_model3 <- rbind(cv_results_model3, tmp_df)
  }
}

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
fig_6c <- ggplot() +
  geom_roc(
    data = cv_results_model1,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = pred
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = "black"
  ) +
  geom_roc(
    data = cv_results_model2,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = pred
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Sunset2", n = 7)[2]
  ) +
  geom_roc(
    data = cv_results_model3,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = pred
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Lake", n = 9)[4]
  ) +
  geom_roc(
    data = cosmis_vs_others,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = - cosmis
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 7)[7]
  ) + 
  geom_roc(
    data = cosmis_vs_others,
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
  filename = paste(data_dir, "fig_6c_6d", "fig_6c.svg", sep = "/"),
  plot = fig_6c,
  width = 5,
  height = 5,
  units = "in",
  device = "svg",
)

# ggplot theme for ROC plot
plot_margin <- margin(
  t = 0.8, r = 0.8, b = 0.5, l = 0.5, unit = "cm"
)
sub_roc_plot_theme <- theme_classic() + theme(
  panel.border = element_rect(
    colour = "black", 
    size = 1, 
    fill = "transparent"
  ),
  axis.text = element_text(
    size = 28, color = "black"
  ),
  axis.title.x = element_text(
    color = "black", size = 32, margin = margin(t = 10, r = 10)
  ),
  axis.title.y = element_text(
    color = "black", size = 32, margin = margin(r = 10)
  ),
  axis.ticks.length = unit(.25, "cm"),
  axis.ticks = element_line(size = 1.0),
  legend.position = "none",
  plot.margin = plot_margin
)


fig_6d <- ggplot() +
  geom_roc(
    data = cv_results_model1,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = pred
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = "black"
  ) +
  geom_roc(
    data = cv_results_model2,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = pred
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Sunset2", n = 7)[2]
  ) +
  geom_roc(
    data = cv_results_model3,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = pred
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Lake", n = 9)[4]
  ) +
  geom_roc(
    data = cosmis_vs_others,
    mapping = aes(
      d = factor(class, levels = c("benign", "pathogenic")),
      m = - cosmis
    ),
    n.cuts = 0,
    labels = FALSE,
    show.legend = FALSE,
    size = 1.2,
    colour = pnw_palette(name = "Shuksan2", n = 7)[7]
  ) + 
  geom_roc(
    data = cosmis_vs_others,
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
    breaks = seq(0, 0.15, 0.05), 
    labels = c("0.0", "0.05", "0.10", "0.15"), 
    limits = c(0.0, 0.15)
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.6, 0.2), 
    labels = c("0.0", "0.20", "0.40", "0.6"), 
    limits = c(0.0, 0.6)
  ) +
  sub_roc_plot_theme

# save the scatter plot to disk
ggsave(
  filename = paste(data_dir, "fig_6c_6d", "fig_6d.svg", sep = "/"),
  plot = fig_6d,
  width = 5,
  height = 8,
  units = "in",
  device = "svg",
)