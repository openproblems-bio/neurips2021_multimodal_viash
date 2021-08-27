cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)

## VIASH START
par <- list(
  input = list.files("src", recursive = TRUE, full.names = TRUE, pattern = "*.tsv"),
  output = "results/meta_metrics.tsv"
)
## VIASH END

metric_paths <- par$input[grepl("/metrics/", par$input)]

cat("Reading input files\n")
metric_meta <- 
  map_df(metric_paths, function(x) {
    x <- gsub("//*", "/", x)
    df <- read_tsv(x, col_types = c(metric_id = "c", metric_min = "f", metric_max = "f", metric_higherisbetter = "l"))

    df$task <- gsub(".*src/([^/]*)/.*", "\\1", x)
    df$component_name <- gsub(".*metrics/([^/]*)/.*", "\\1", x)

    df
  })


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

cat("Writing output file\n")
write_tsv(metric_meta %>% select(task, component_name, everything()), par$output)
