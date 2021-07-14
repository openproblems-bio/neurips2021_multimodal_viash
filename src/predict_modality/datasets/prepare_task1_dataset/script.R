cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input_mod1 = "resources_test/common/pbmc_1k_protein_v3.normalize.output_rna.h5ad",
  input_mod2 = "resources_test/common/pbmc_1k_protein_v3.normalize.output_mod2.h5ad",
  output_mod1 = "output_mod1.h5ad",
  output_mod2 = "output_mod2.h5ad",
  output_solution = "solution.h5ad",
  seed = 1L,
  max_mod2_columns = 5000L
)
## VIASH END

cat("Reading input data\n")
ad1_raw <- anndata::read_h5ad(par$input_mod1)
ad2_raw <- anndata::read_h5ad(par$input_mod2)
ad1_mod <- unique(ad1_raw$var[["feature_types"]])
ad2_mod <- unique(ad2_raw$var[["feature_types"]])
new_dataset_id <- paste0(ad1_raw$uns[["dataset_id"]], "_task1_", tolower(ad1_mod))
common_uns <- list(dataset_id = new_dataset_id)

cat("Determining train/test split\n")
split <- 
  if (!is.null(ad1_raw$obs[["experiment"]]) && all(ad1_raw$obs[["experiment"]] %in% c("train", "test"))) {
    ad1_raw$obs[["experiment"]]
  } else {
    set.seed(par$seed)
    ix <- sample.int(
      nrow(ad1_raw), 
      size = nrow(ad1_raw) * 0.66,
      replace = FALSE
    )
    ifelse(seq_len(nrow(ad1_raw)) %in% ix, "train", "test")
  }
splor <- order(split)
ad1_X <- ad1_raw$X[splor,, drop = FALSE]
ad2_X <- ad2_raw$X[splor,, drop = FALSE]
split <- split[splor]

if (!is.null(par$max_mod2_columns) && par$max_mod2_columns < ncol(ad2_X)) {
  cat("Sampling mod2 columns\n")
  ad2_ix <- sample.int(ncol(ad2_X), par$max_mod2_columns)
  ad2_X <- ad2_X[, ad2_ix, drop = FALSE]
}

cat("Creating mod1 object\n")
out_mod1 <- anndata::AnnData(
  X = ad1_X,
  var = data.frame(
    row.names = colnames(ad1_X),
    feature_types = rep(ad1_mod, ncol(ad1_X))
  ),
  obs = data.frame(
    row.names = rownames(ad1_X),
    split = split
  ),
  uns = common_uns
)

cat("Creating mod2 object\n")
ad2_var <- data.frame(
  row.names = colnames(ad2_X),
  feature_types = rep(ad2_mod, ncol(ad2_X))
)
ad2_obs <- data.frame(
  row.names = rownames(ad2_X),
  split = split
)

out_mod2 <- anndata::AnnData(
  X = ad2_X[split == "train",, drop = FALSE],
  var = ad2_var,
  obs = ad2_obs[split == "train",, drop = FALSE],
  uns = common_uns
)

cat("Create solution object\n")
out_solution <- anndata::AnnData(
  X = ad2_X[split == "test",, drop = FALSE],
  var = ad2_var,
  obs = ad2_obs[split == "test",, drop = FALSE],
  uns = common_uns
)

cat("Saving output files as h5ad\n")
zzz <- out_mod1$write_h5ad(par$output_mod1, compression = "gzip")
zzz <- out_mod2$write_h5ad(par$output_mod2, compression = "gzip")
zzz <- out_solution$write_h5ad(par$output_solution, compression = "gzip")