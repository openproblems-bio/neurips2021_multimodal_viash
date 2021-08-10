cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(GA, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input_rna = "resources_test/common/test_resource.output_rna.h5ad",
  input_mod2 = "resources_test/common/test_resource.output_mod2.h5ad",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad"
)
## VIASH END

ideal_train_pct <- 0.66

cat("Reading h5ad file\n")
ad_rna <- anndata::read_h5ad(par$input_rna)
ad_mod2 <- anndata::read_h5ad(par$input_mod2)

cat("Determining train/test group\n")
is_train <-
  if (!is.null(ad_rna$obs[["batch"]])) {
    # use batch label to optimise the best train/test split
    batch <- as.character(ad_rna$obs[["batch"]])
    nums <- table(batch)
    nums <- setNames(as.numeric(nums), names(nums))

    # try to find a split that best matches 0.66
    ga_out <- GA::ga(
      type = "binary",
      nBits = length(nums),
      fitness = function(x) {
        -abs(sum(nums[x == 1]) / sum(nums) - ideal_train_pct)
      },
      maxiter = 500,
      monitor = FALSE
    )
    train_batches <- names(nums)[ga_out@solution[1,] == 1]
    batch %in% train_batches
  } else {
    # TODO: should batch effect be simulated
    # if there is no batch in the dataset??
    
    # split randomly
    set.seed(par$seed)
    ix <- sample.int(
      nrow(ad_rna), 
      size = nrow(ad_rna) * ideal_train_pct,
      replace = FALSE
    )
    seq_len(nrow(ad_rna)) %in% ix
  }

ad_rna.obs["is_train"] <- is_train
ad_mod2.obs["is_train"] <- is_train

cat("Writing mod1 data\n")
print(ad_rna)
zzz <- ad_rna$write_h5ad(par$output_rna, compression = "gzip")

cat("Writing mod2 data\n")
print(ad_mod2)
zzz <- ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")
