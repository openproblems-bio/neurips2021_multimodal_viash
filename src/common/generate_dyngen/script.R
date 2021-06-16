## VIASH START
par <- list(
  backbone = "bifurcating",
  num_cells = 100,
  num_genes = 120,
  num_simulations = 3,
  num_threads = 1,
  model = NULL,
  output = "output_dyngen.h5ad",
  output_censored = "output_dyngen_censored.h5ad",
  plot = "output_dyngen_plot.pdf",
  ssa_tau = 30 / 3600,
  census_interval = 1,
  compute_cellwise_grn = FALSE,
  compute_rna_velocity = FALSE
)
## VIASH END

options(tidyverse.quiet = TRUE)
library(tidyverse)
library(dyngen, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

# determine backbone
backbones <- list_backbones()

backbone <- backbones[[par$backbone]]()

# generate initial config
num_tfs <- nrow(backbone$module_info)
num_targets <- ceiling(0.8 * (par$num_genes - num_tfs))
num_hks <- par$num_genes - num_tfs - num_targets

init_model <- initialise_model(
  backbone = backbone,
  num_cells = par$num_cells,
  num_tfs = num_tfs,
  num_targets = num_targets,
  num_hks = num_hks,
  simulation_params = simulation_default(
    census_interval = par$census_interval,
    ssa_algorithm = ssa_etl(tau = par$ssa_tau),
    experiment_params = simulation_type_wild_type(
      num_simulations = par$num_simulations
    ),
    compute_cellwise_grn = par$compute_cellwise_grn,
    compute_rna_velocity = par$compute_rna_velocity
  ),
  num_cores = par$num_threads,
  verbose = TRUE
)

# run simulations
out <- generate_dataset(
  model = init_model,
  format = "anndata",
  make_plots = TRUE
)

if (!is.null(par$model)) {
  write_rds(out$model, par$model, compress = "gz")
}

dataset <- out$dataset

dataset$uns[["dataset_id"]] <- par$id

# if a censored version of the dataset needs to be generated
if (!is.null(par$output_censored)) {
  counts <- dataset$X
  counts_protein <- dataset$layers["counts_protein"]

  # shuffle
  shuffle_cells <- sample.int(nrow(counts_protein))
  shuffle_genes <- sample.int(ncol(counts_protein))

  counts <- counts[, shuffle_genes, drop = FALSE]
  counts_protein <- counts_protein[shuffle_cells, shuffle_genes, drop = FALSE]

  # rename genes, remove cell names
  rownames(counts) <- rownames(counts_protein) <- NULL
  colnames(counts) <- colnames(counts_protein) <- paste0("gene_", seq_len(ncol(counts)))

  # create new dataset object
  dataset_censored <- anndata::AnnData(
    X = counts,
    layers = list(
      protein = counts_protein
    ),
    uns = list(
      dataset_id = par$dataset_id
    )
  )

  dataset$uns$shuffle_cells <- shuffle_cells-1
  dataset$uns$shuffle_genes <- shuffle_genes-1

  dataset_censored$write_h5ad(par$output_censored, compression = "gzip")
}


dataset$write_h5ad(par$output, compression = "gzip")

# save plot (if need be)
if (!is.null(par$plot)) {
  ggsave(par$plot, out$plot, width = 20, height = 16)
}
