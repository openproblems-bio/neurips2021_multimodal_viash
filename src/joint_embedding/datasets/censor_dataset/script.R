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
  output_solution = "solution.h5ad"
)
## VIASH END

cat("Reading input data\n")
ad1_raw <- anndata::read_h5ad(par$input_mod1)
ad2_raw <- anndata::read_h5ad(par$input_mod2)
common_uns <- list(
  dataset_id = paste0(ad1_raw$uns[["dataset_id"]], "_JE"),
  organism = ad1_raw$uns[["organism"]]
)

cat("Creating mod1 object\n")
out_mod1 <- anndata::AnnData(
  X = ad1_raw$X,
  var = ad1_raw$var %>% select(one_of("gene_ids"), feature_types),
  obs = ad1_raw$obs %>% select(one_of("batch")),
  uns = common_uns
)

cat("Creating mod2 object\n")
out_mod2 <- anndata::AnnData(
  X = ad2_raw$X,
  var = ad2_raw$var %>% select(one_of("gene_ids"), feature_types),
  obs = ad2_raw$obs %>% select(one_of("batch")),
  uns = common_uns
)

cat("Create solution object\n")
out_solution <- anndata::AnnData(
  X = ad1_raw$X,
  var = ad1_raw$var %>% select(one_of("gene_ids"), feature_types),
  obs = ad1_raw$obs %>% select(
    one_of("batch", "cell_type", "pseudotime_order_GEX", "pseudotime_order_ATAC", "pseudotime_order_ADT", "S_score", "G2M_score")
  ),
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
