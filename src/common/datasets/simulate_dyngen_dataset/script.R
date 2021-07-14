## VIASH START
par <- list(
  id = "test",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad",
  backbone = "linear",
  num_cells = 100,
  num_genes = 120,
  num_simulations = 10,
  num_threads = 10,
  plot = "plot.pdf",
  ssa_tau = 30 / 3600,
  census_interval = 1,
  store_chromatin = FALSE,
  store_protein = TRUE,
  num_proteins = 50
)
## VIASH END

if (par$store_protein == par$store_chromatin) {
  cat("Warning: Strictly pass one of --store_protein and --store_chromatin, not neither or both.\n")
}

cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(dyngen, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

cat("Creating dyngen backbone\n")
backbones <- list_backbones()

if (is.null(par$backbone)) {
  par$backbone <- sample(names(backbones), 1)
}

backbone <- backbones[[par$backbone]]()

cat("Generating regulatory network\n")
num_tfs <- nrow(backbone$module_info)
num_targets <- ceiling(0.8 * (par$num_genes - num_tfs))
num_hks <- par$num_genes - num_tfs - num_targets
num_cells_train <- round(par$num_cells * .66)
num_cells_test <- par$num_cells - num_cells_train

model_init <- initialise_model(
  backbone = backbone,
  num_cells = num_cells_train,
  num_tfs = num_tfs,
  num_targets = num_targets,
  num_hks = num_hks,
  simulation_params = simulation_default(
    census_interval = par$census_interval,
    ssa_algorithm = ssa_etl(tau = par$ssa_tau),
    experiment_params = simulation_type_wild_type(
      num_simulations = round(0.66*par$num_simulations)
    ),
    compute_cellwise_grn = par$store_chromatin
  ),
  num_cores = par$num_threads,
  verbose = FALSE
) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics()

cat("Running simulations for training cells\n")
model_train <-
  model_init %>%
  generate_gold_standard() %>%
  generate_cells() %>%
  generate_experiment()

cat("Running simulations for test cells\n")
model_init$num_cells <-
  model_init$numbers$num_cells <-
  num_cells_test
model_init$simulation_params$experiment_params <-
  simulation_type_wild_type(round(0.37*par$num_simulations))
model_test <-
  model_init %>%
  generate_gold_standard() %>%
  generate_cells() %>%
  generate_experiment()

cat("Combine simulations into one dataset\n")
model <- combine_models(
  list(train = model_train, test = model_test),
  duplicate_gold_standard = FALSE
)
dataset <- as_list(model)

cat("Create RNA dataset\n")
obs <- dataset$cell_info %>% 
  rename(experiment = model) %>%
  column_to_rownames("cell_id")
var <- dataset$feature_info %>% 
  select(feature_id, module_id, basal, burn, independence, color, is_tf, is_hk) %>% 
  column_to_rownames("feature_id")
ad_mod1 <- anndata::AnnData(
  X = dataset$counts,
  obs = obs,
  var = var %>% mutate(feature_types = "GEX"),
  uns = list(
    dataset_id = par$id
  )
)

if (par$store_protein) {
  cat("Processing Antibody data\n")
  # construct AbSeq-like data from protein counts
  # TODO: use real AbSeq data to map distributions
  counts_protein <- dataset$counts_protein
  var_protein <- var %>% mutate(feature_types = "ADT")

  # sample 50 genes
  if (ncol(counts_protein) > par$num_proteins) {
    sample_genes <- sample.int(ncol(counts_protein), par$num_proteins)
    counts_protein <- counts_protein[,sample_genes, , drop = FALSE]
    var_protein <- var_protein[sample_genes, , drop = FALSE]
  }

  ad_mod2 <- anndata::AnnData(
    X = counts_protein,
    obs = obs,
    var = var_protein,
    uns = list(
      dataset_id = par$id
    )
  )
}

if (par$store_chromatin) {
  cat("Processing ATAC data\n")
  # constuct atac-like data from single cell regulatory network
  # TODO: use real atac data to map distributions
  mat <- dataset$regulatory_network_sc %>%
    mutate(
      edge = factor(paste0(as.character(regulator), "->", as.character(target)))
    )
  regsc <- Matrix::sparseMatrix(
    i = as.integer(mat$cell_id),
    j = as.integer(mat$edge),
    x = pmax(mat$strength, 0)*100
  )
  rownames(regsc) <- dataset$cell_ids
  colnames(regsc) <- paste0("region_", seq_len(ncol(regsc)))
  var_atac <- data.frame(
    row.names = colnames(regsc),
    feature_types = rep("ATAC", ncol(regsc))
  )

  ad_mod2 <- anndata::AnnData(
    X = regsc,
    obs = obs,
    var = var_atac,
    uns = list(
      dataset_id = par$id
    )
  )
}

cat("Write h5ad files\n")
ad_mod1$write_h5ad(par$output_rna, compression = "gzip")
ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")

if (!is.null(par$plot)) {
  cat("Storing summary plot\n")
  g <- plot_summary(model)
  ggsave(par$plot, g, width = 20, height = 16)
}



