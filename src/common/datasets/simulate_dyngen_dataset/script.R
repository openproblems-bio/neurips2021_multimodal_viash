cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(dyngen, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  id = "test",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad",
  num_cells = 300,
  num_genes = 120,
  num_simulations = 50,
  num_threads = 10,
  plot = "plot.pdf",
  ssa_tau = 30 / 3600,
  census_interval = 1,
  store_chromatin = FALSE,
  store_protein = TRUE,
  num_proteins = 50,
  cache_dir = tools::R_user_dir("dyngen", "data"),
  seed = 1
)
## VIASH END

if (!is.null(par$seed)) {
  set.seed(par$seed)
}
if (par$store_protein == par$store_chromatin) {
  cat("Warning: Strictly pass one of --store_protein and --store_chromatin, not neither or both.\n")
}

# start from given backbone and a cyclic backbone
backbone_init <- bblego(
  bblego_start("A", type = "simple"),
  bblego_linear("A", "B", num_modules = 10),
  bblego_linear("B", "C", num_modules = 10),
  bblego_linear("C", "D", num_modules = 10),
  bblego_end("D")
)


module_info <- bind_rows(
  backbone_init$module_info %>% select(-color),
  tribble(
    ~module_id, ~basal, ~burn, ~independence,
    "CC1", 0, TRUE, 1,
    "CC2", 1, TRUE, 1,
    "CC3", 0, TRUE, 1,
    "CC4", 1, TRUE, 1,
    "CC5", 0, TRUE, 1
  )
)

module_network <- bind_rows(
  backbone_init$module_network,
  tribble(
    ~from, ~to, ~effect, ~strength, ~hill,
    "Burn1", "CC1", 1L, 1, 2,
    "CC1", "CC2", -1L, 100, 2,
    "CC2", "CC3", 1L, 1, 2,
    "CC3", "CC4", -1L, 100, 2,
    "CC4", "CC5", 1L, 1, 2,
    "CC5", "CC1", -1L, 100, 2
  )
)

backbone <- backbone(
  module_info = module_info,
  module_network = module_network,
  expression_patterns = backbone_init$expression_patterns
)

burn_time <- simtime_from_backbone(backbone, burn = TRUE) * 4
total_time <- simtime_from_backbone(backbone, burn = FALSE) / 3

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
    burn_time = burn_time,
    total_time = total_time,
    census_interval = par$census_interval,
    ssa_algorithm = ssa_etl(tau = par$ssa_tau),
    experiment_params = simulation_type_wild_type(
      num_simulations = round(0.66*par$num_simulations)
    ),
    compute_cellwise_grn = par$store_chromatin
  ),
  num_cores = par$num_threads,
  verbose = TRUE,
  download_cache_dir = par$cache_dir
) %>%
  generate_tf_network() %>%
  generate_feature_network()

cat("Running simulations for training cells\n")
model_train <-
  model_init %>%
  generate_kinetics() %>%
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
  generate_kinetics() %>% # use different kinetics
  generate_gold_standard() %>%
  generate_cells() %>%
  generate_experiment()

cat("Combine simulations into one dataset\n")
model <- combine_models(
  list(
    batch1 = model_train,
    batch2 = model_test
  ),
  duplicate_gold_standard = FALSE
)
dataset <- as_dyno(model) %>% 
  dynwrap::simplify_trajectory()

milestone_labels <- c("sA" = "start", "sEndD" = "end")
# check whether output dataset looks nice
# and whether batch effect is present
# plot_summary(model)


cat("Compute cell cycle scores\n")
phase_scores <- map_df(unique(dataset$cell_info$model), function(mod) {
  ix <- dataset$cell_info$model == mod
  cp <- dataset$counts_protein[ix, ]
  expr1 <- dynutils::scale_quantile(cp[, "CC2_TF1"])
  expr2 <- dynutils::scale_quantile(cp[, "CC5_TF1"])
  tibble(
    cell_id = rownames(cp),
    S_score = (expr1 - expr2) / 2 + .5,
    G2M_score = (expr2 - expr1) / 2 + .5
  )
})

# ggplot(phase_scores) + geom_point(aes(S_score, G2M_score))
cat("Process trajectory percentages")
trajectory_topology <-
  dataset$milestone_network %>%
  transmute(
    from = milestone_labels[from],
    to = milestone_labels[to],
    length
  )
trajectory_percentages <- 
  dataset$milestone_percentages %>%
  mutate(
    milestone_id = milestone_labels[milestone_id]
  ) %>%
  reshape2::acast(cell_id~milestone_id, value.var = "percentage") %>%
  Matrix::Matrix(sparse = TRUE)

trajectory_percentages <- trajectory_percentages[dataset$cell_ids, unname(milestone_labels)]

cat("Create RNA dataset\n")
celltypes <- dataset$milestone_percentages %>%
  group_by(cell_id) %>%
  slice(which.max(percentage)) %>%
  ungroup() %>%
  select(cell_id, cell_type = milestone_id)
obs <- dataset$cell_info %>%
  left_join(celltypes, by = "cell_id") %>%
  left_join(phase_scores, by = "cell_id") %>%
  rename(batch = model) %>%
  column_to_rownames("cell_id")

# ggplot(obs) + geom_point(aes(pseudotime, G2M_score, colour = batch))
var <- dataset$feature_info %>%
  select(feature_id, module_id, basal, burn, independence, color, is_tf, is_hk) %>%
  column_to_rownames("feature_id")
ad_mod1 <- anndata::AnnData(
  X = dataset$counts,
  obs = obs,
  var = var %>% mutate(feature_types = "GEX"),
  obsm = list(
    trajectory_percentages = trajectory_percentages
  ),
  uns = list(
    dataset_id = par$id,
    trajectory_topology = trajectory_topology
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
    counts_protein <- counts_protein[, sample_genes, , drop = FALSE]
    var_protein <- var_protein[sample_genes, , drop = FALSE]
  }

  ad_mod2 <- anndata::AnnData(
    X = counts_protein,
    obs = obs,
    var = var_protein,
    obsm = list(
      trajectory_percentages = trajectory_percentages
    ),
    uns = list(
      dataset_id = par$id,
      trajectory_topology = trajectory_topology
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
    x = pmax(mat$strength, 0) * 100
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
    obsm = list(
      trajectory_percentages = trajectory_percentages
    ),
    uns = list(
      dataset_id = par$id,
      organism = "synthetic"
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
