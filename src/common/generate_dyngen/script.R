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
  compute_atac = TRUE,
  compute_rna_velocity = FALSE,
  store_protein = TRUE
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
    compute_cellwise_grn = par$compute_atac,
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

# get mrna counts
counts <- out$dataset$X

# construct Ab-like data from protein counts
# TODO: use real Ab data to map distributions
counts_ab <- out$dataset$layers[["counts_protein"]]

# sample 50 genes
sample_genes <-
  if (ncol(counts_ab) > 50) {
    sample.int(ncol(counts_ab), 50)
  } else {
    seq_len(ncol(counts_ab))
  }

counts_ab[,-sample_genes] <- 0

counts_ab <- Matrix::drop0(counts_ab)

# constuct atac-like data from single cell regulatory network
# TODO: use real atac data to map distributions
regulatory_network <- out$dataset$uns[["regulatory_network"]]
regulatory_network_sc <- out$dataset$obsm[["regulatory_network_sc"]]

library(Matrix, quietly = TRUE)
summ <- summary(regulatory_network_sc)
regnet_ix <- tibble(
  cell_i = summ$i,
  reg_i = summ$j,
  reg_value = summ$x
) %>%
  left_join(
    regulatory_network %>%
      transmute(
        reg_i = row_number(),
        gene_i = match(target, colnames(counts))
      ),
    by = "reg_i"
  )

atac_ix <-
  regnet_ix %>%
  group_by(cell_i, gene_i) %>%
  summarise(
    atac = sum(reg_value),
    .groups = "drop"
  )
atac <- Matrix::sparseMatrix(
  i = atac_ix$cell_i,
  j = atac_ix$gene_i,
  x = pmax(atac_ix$atac, 0),
  dims = dim(counts),
  dimnames = dimnames(counts)
)

adata <- anndata::AnnData(
  X = NULL,
  obs = out$dataset$obs,
  var = out$dataset$var,
  shape = dim(counts),
  layers = list(
    mrna = counts,
    antibody = counts_ab,
    atac = atac
  ),
  uns = list(
    dataset_id = par$dataset_id
  )
)

adata$write_h5ad(par$output, compression = "gzip")

# save plot (if need be)
if (!is.null(par$plot)) {
  ggsave(par$plot, out$plot, width = 20, height = 16)
}



