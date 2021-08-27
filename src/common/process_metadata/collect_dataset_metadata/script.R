cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)

## VIASH START
par <- list(
  input_tsvs = list.files("src/common/datasets/", recursive = TRUE, full.names = TRUE, pattern = "*.tsv"),
  input_h5ads = list.files("output/public_datasets/common/", recursive = TRUE, full.names = TRUE, pattern = "*.h5ad"),
  output = "results/meta_datasets.tsv"
)
## VIASH END

cat("Reading input files\n")
dataset_metadata_tsv <- 
  map_df(par$input_tsvs, read_tsv) %>% 
  select(dataset_id = id, geo_id)

dataset_metadata_h5ad <- 
  map_df(par$input_h5ads, function(x) {
    ad <- anndata::read_h5ad(x, backed = TRUE)
    tibble(
      dataset_id = ad$uns$dataset_id,
      organism = ad$uns$organism,
      num_cells = nrow(ad),
      num_features = ncol(ad),
      modality = unique(ad$var$feature_types),
      num_batches = length(unique(ad$obs$batch)),
      num_cell_types = length(unique(ad$obs$cell_type))
    )
  })

combined <- 
  dataset_metadata_tsv %>% 
  left_join(dataset_metadata_h5ad %>% filter(modality == "GEX") %>% rename(num_features_gex = num_features), by = "dataset_id") %>%
  left_join(dataset_metadata_h5ad %>% filter(modality != "GEX") %>% select(dataset_id, num_features_mod2 = num_features), by = "dataset_id") %>%
  filter(!is.na(num_cells))

cat("Writing output file\n")
write_tsv(combined, par$output)
