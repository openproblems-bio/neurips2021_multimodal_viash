## VIASH START
par <- list(
  backbone = "bifurcating",
  num_cells = 100,
  num_genes = 120,
  num_simulations = 3,
  num_threads = 3,
  output = "output.h5ad",
  plot = "plot.pdf",
  ssa_tau = 30 / 3600,
  census_interval = 1,
  store_chromatin = TRUE,
  store_rna_velocity = FALSE,
  store_protein = TRUE,
  num_proteins = 50
)
## VIASH END

if (par$store_protein == par$store_chromatin) {
  cat("Warning: Strictly pass one of --store_protein and --store_chromatin, not neither or both.\n")
}

options(tidyverse.quiet = TRUE)
library(tidyverse)
library(dyngen, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

# determine backbone
backbones <- list_backbones()

if (is.null(par$backbone)) {
  par$backbone <- sample(names(backbones), 1)
}

backbone <- backbones[[par$backbone]]()

# generate initial config
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
      num_simulations = par$num_simulations
    ),
    compute_cellwise_grn = par$store_chromatin,
    compute_rna_velocity = par$store_rna_velocity
  ),
  num_cores = par$num_threads,
  verbose = TRUE
) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics()

# run train simulations
model_train <-
  model_init %>%
  generate_gold_standard() %>%
  generate_cells() %>%
  generate_experiment()

# run test simulations
model_init$num_cells <-
  model_init$numbers$num_cells <-
  num_cells_test
model_test <-
  model_init %>%
  generate_gold_standard() %>%
  generate_cells() %>%
  generate_experiment()

# combine into one dataset
model <- combine_models(
  list(train = model_train, test = model_test),
  duplicate_gold_standard = FALSE
)
dataset <- as_anndata(model)


# create adata dataset
counts <- dataset$X

adata <- anndata::AnnData(
  X = counts,
  obs = dataset$obs %>% rename(experiment = model),
  var = dataset$var %>% select(module_id, basal, burn, independence, color, is_tf, is_hk),
  uns = list(
    dataset_id = par$id
  )
)

if (par$store_protein) {
  # construct AbSeq-like data from protein counts
  # TODO: use real AbSeq data to map distributions
  counts_protein <- dataset$layers[["counts_protein"]]

  # sample 50 genes
  if (ncol(counts_protein) > par$num_proteins) {
    sample_genes <- sample.int(ncol(counts_protein), par$num_proteins)
    counts_protein <- counts_protein[,sample_genes]
  }

  adata$obsm[["protein"]] <- counts_protein
  adata$uns[["protein_varnames"]] <- colnames(counts_protein)
}

if (par$store_chromatin) {
  # constuct atac-like data from single cell regulatory network
  # TODO: use real atac data to map distributions
  regulatory_network <- dataset$uns[["regulatory_network"]]
  regulatory_network_sc <- dataset$obsm[["regulatory_network_sc"]]

  adata$obsm[["chromatin"]] <- regulatory_network_sc
  adata$uns[["chromatin_colnames"]] <- paste0("region_", seq_len(ncol(regulatory_network_sc)))
}

adata$write_h5ad(par$output, compression = "gzip")

# save plot (if need be)
if (!is.null(par$plot)) {
  g <- plot_summary(model)
  ggsave(par$plot, g, width = 20, height = 16)
}



