cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)

## VIASH START
par <- list(
  input = c(
  "resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.scores.h5ad", 
  "resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.scores.h5ad"
  ),
  output = "tmp/task1_scores.tsv"
)
## VIASH END

cat("Reading input files\n")
output <- map_df(par$input, read_tsv)

cat("Writing output file\n")
write_tsv(output, par$output)