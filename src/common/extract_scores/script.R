cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  # input = "resources_test/predict_modality/test_resource.scores.h5ad",
  input = list.files("work/f2/76020c558a40e0d469f18dc56bb258", pattern = "*.(h5ad|output)$", full.names = TRUE),
  output = "output/pilot/joint_embedding/output.extract_scores.output.tsv",
  summary = "output/pilot/joint_embedding/output.extract_scores.summary.tsv",
  method_meta = NULL,
  dataset_meta = "work/f2/76020c558a40e0d469f18dc56bb258/solutions_meta.tsv",
  metric_meta = "work/f2/76020c558a40e0d469f18dc56bb258/meta.bind_tsv_rows.output.tsv"
)
par$input <- par$input[!duplicated(basename(par$input))]
inp <- par$input[[1]]
## VIASH END

cat("Reading score meta files\n")
metric_meta <- map_df(
  par$metric_meta,
  read_tsv,
  col_types = cols(
    metric_id = "c",
    metric_min = "d",
    metric_max = "d",
    metric_higherisbetter = "l"
  )
)
# prepend 'metric_' to colnames if this was not done already
colnames(metric_meta) <- paste0("metric_", gsub("^metric_", "", colnames(metric_meta)))

cat("Reading input h5ad files\n")
scores <- map_df(par$input, function(inp) {
  cat("Reading '", inp, "'\n", sep = "")
  ad <- anndata::read_h5ad(inp)

  for (uns_name in c("dataset_id", "method_id", "metric_ids", "metric_values")) {
    expect_true(
      uns_name %in% names(ad$uns),
      info = paste0("File ", inp, " must contain `uns['", uns_name, "']`")
    )
  }

  as_tibble(ad$uns[c("dataset_id", "method_id", "metric_ids", "metric_values")])
}) %>%
  rename(metric_id = metric_ids, value = metric_values)

expect_true(
  all(unique(scores$metric_id) %in% metric_meta$metric_id),
  info = "All metric_ids in h5ad should be mentioned in the metric_meta"
)

# create method meta
method_meta <-
  if (is.null(par$method_meta)) {
    cat("No method meta found, defaulting to all method ids present in the h5ads.\n")
    tibble(
      method_id = unique(scores$method_id)
    )
  } else {
    cat("Reading method meta files\n")
    map_df(
      par$method_meta,
      read_tsv,
      col_types = cols(
        method_id = "c",
        .default = "?"
      )
    )
  }
print(method_meta)
# prepend 'method_' to colnames if this was not done already
colnames(method_meta) <- paste0("method_", gsub("^method_", "", colnames(method_meta)))

# create method meta
dataset_meta <- 
  if (is.null(par$dataset_meta)) {
    cat("No dataset meta found, defaulting to all dataset ids present in the h5ads.\n")
    tibble(
      dataset_id = unique(scores$dataset_id)
    )
  } else {
    cat("Reading method meta files\n")
    map_df(
      par$dataset_meta,
      read_tsv,
      col_types = cols(
        dataset_id = "c",
        .default = "?"
      )
    )
  }
# prepend 'dataset_' to colnames if this was not done already
colnames(dataset_meta) <- paste0("dataset_", gsub("^dataset_", "", colnames(dataset_meta)))

cat("Creating default scores for missing entries based on metrics meta\n")
missing_values <- crossing(
  dataset_meta %>% select(dataset_id),
  method_meta %>% select(method_id),
  metric_meta %>% transmute(
    metric_id = metric_id,
    value = ifelse(metric_higherisbetter, metric_min, metric_max)
  )
)

cat("Creating final scores object\n")
final_scores <- bind_rows(
  scores %>%
    inner_join(method_meta, by = "method_id") %>%
    inner_join(dataset_meta, by = "dataset_id") %>%
    inner_join(metric_meta %>% select(metric_id), by = "metric_id"),
  anti_join(missing_values, scores, by = c("method_id", "dataset_id", "metric_id"))
)
cat("Writing results to tsv files\n")
write_tsv(final_scores, par$output)
