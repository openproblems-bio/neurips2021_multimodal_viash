cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)

## VIASH START
par <- list(
  input = list.files("output/datasets/predict_modality/", recursive = TRUE, full.names = TRUE, pattern = "output_test_mod2.h5ad"),
  output = "results/meta_datasets.tsv"
)
## VIASH END

cat("Reading input files\n")
output <-
  map_df(par$input, function(x) {
    cat("Reading ", x, "\n", sep = "")
    ad <- anndata::read_h5ad(x)
    df <- tibble(
      dataset_id = ad$uns$dataset_id,
      modality = unique(ad$var$feature_types),
      default_rmse = sqrt(sum(ad$X@x^2) / length(ad$X)),
      default_mae = sum(abs(ad$X@x)) / length(ad$X)
    )
    rm(ad)
    df
  })

cat("Writing output file\n")
write_tsv(output, par$output)
