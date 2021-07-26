cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input_mod1 = "resources_test/common/test_resource.output_rna.h5ad",
  input_mod2 = "resources_test/common/test_resource.output_mod2.h5ad",
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
new_dataset_id <- paste0(ad1_raw$uns[["dataset_id"]], "_predmod_", tolower(ad1_mod), "2", tolower(ad2_mod))
common_uns <- list(dataset_id = new_dataset_id)

cat("Determining train/test group\n")
group <- 
  if (!is.null(ad1_raw$obs[["batch"]]) && all(ad1_raw$obs[["batch"]] %in% c("train", "test"))) {
    ad1_raw$obs[["batch"]]
  } else if (!is.null(ad1_raw$obs[["batch"]])) {
    # use batch label to optimise the best train/test split
    batch <- as.character(ad1_raw$obs[["batch"]])
    nums <- table(batch)
    nums <- setNames(as.numeric(nums), names(nums))

    # try to find a split that best matches 0.66
    ga_out <- GA::ga(
      type = "binary",
      nBits = length(nums),
      fitness = function(x) {
        -abs(sum(nums[x == 1]) / sum(nums) - 0.66)
      },
      maxiter = 500,
      monitor = FALSE
    )
    train_batches <- names(nums)[ga_out@solution[1,] == 1]
    ifelse(batch %in% train_batches, "train", "test")
  } else {
    # split randomly
    set.seed(par$seed)
    ix <- sample.int(
      nrow(ad1_raw), 
      size = nrow(ad1_raw) * 0.66,
      replace = FALSE
    )
    ifelse(seq_len(nrow(ad1_raw)) %in% ix, "train", "test")
  }
splor <- order(group)
ad1_X <- ad1_raw$X[splor,, drop = FALSE]
ad2_X <- ad2_raw$X[splor,, drop = FALSE]
group <- group[splor]

if (!is.null(par$max_mod1_columns) && par$max_mod1_columns < ncol(ad1_X)) {
  cat("Sampling mod1 columns\n")
  ad1_ix <- sample.int(ncol(ad1_X), par$max_mod1_columns)
  ad1_X <- ad1_X[, ad1_ix, drop = FALSE]
}
if (!is.null(par$max_mod2_columns) && par$max_mod2_columns < ncol(ad2_X)) {
  cat("Sampling mod2 columns\n")
  ad2_ix <- sample.int(ncol(ad2_X), par$max_mod2_columns)
  ad2_X <- ad2_X[, ad2_ix, drop = FALSE]
}

cat("Creating mod1 object\n")
out_mod1 <- anndata::AnnData(
  X = ad1_X,
  var = ad1_raw$var %>% select(one_of("gene_ids"), feature_types),
  obs = ad1_raw$obs[splor, ] %>% select(one_of("batch")) %>% mutate(group),
  uns = common_uns
)

cat("Creating mod2 object\n")
ad2_var <- ad2_raw$var %>% select(one_of("gene_ids"), feature_types)
ad2_obs <- ad2_raw$obs[splor, ] %>% select(one_of("batch")) %>% mutate(group)

out_mod2 <- anndata::AnnData(
  X = ad2_X[group == "train",, drop = FALSE],
  var = ad2_var,
  obs = ad2_obs[group == "train",, drop = FALSE],
  uns = common_uns
)

cat("Create solution object\n")
out_solution <- anndata::AnnData(
  X = ad2_X[group == "test",, drop = FALSE],
  var = ad2_var,
  obs = ad2_obs[group == "test",, drop = FALSE],
  uns = common_uns
)

cat("Saving output files as h5ad\n")
cat("output_mod1:")
print(out_mod1)
zzz <- out_mod1$write_h5ad(par$output_mod1, compression = "gzip")

cat("output_mod2:")
print(out_mod2)
zzz <- out_mod2$write_h5ad(par$output_mod2, compression = "gzip")

cat("output_solution:")
print(out_solution)
zzz <- out_solution$write_h5ad(par$output_solution, compression = "gzip")