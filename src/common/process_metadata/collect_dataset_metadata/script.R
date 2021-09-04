cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)

## VIASH START
par <- list(
  input_tsvs = list.files("src/common/datasets/", recursive = TRUE, full.names = TRUE, pattern = "*.tsv"),
  input_h5ads = list.files("output/datasets/common/", recursive = TRUE, full.names = TRUE, pattern = "*.h5ad"),
  output = "results/meta_datasets.tsv"
)
## VIASH END

cat("Reading input files\n")
dataset_metadata_tsv <- 
  map_df(par$input_tsvs, read_tsv) %>% 
  select(dataset_id = id, geo_id)

dataset_metadata_h5ad <- 
  map_df(par$input_h5ads, function(x) {
    cat("Reading ", x, "\n", sep = "")
    ad <- anndata::read_h5ad(x, backed = TRUE)
    num_batches <-
      if (ad$uns[["organism"]] == "synthetic" || ad$uns[["batch_type"]] == "real") {
        length(unique(ad$obs$batch))
      } else {
        1
      }
    tibble(
      dataset_id = ad$uns$dataset_id,
      organism = ad$uns$organism,
      num_cells = nrow(ad),
      num_features = ncol(ad),
      modality = unique(ad$var$feature_types),
      num_batches = num_batches,
      num_cell_types = length(unique(ad$obs$cell_type))
    )
  })

combined <-
  dataset_metadata_tsv %>% 
  full_join(dataset_metadata_h5ad %>% filter(modality == "GEX") %>% rename(num_features_mod1 = num_features, mod1_modality = modality), by = "dataset_id") %>%
  full_join(dataset_metadata_h5ad %>% filter(modality != "GEX") %>% select(dataset_id, num_features_mod2 = num_features, mod2_modality = modality), by = "dataset_id") %>%
  filter(!is.na(num_cells)) %>%
  mutate(source = gsub("_.*", "", dataset_id)) %>%
  select(source, dataset_id, geo_id, everything())

cat("Writing output file\n")
write_tsv(combined, par$output)
