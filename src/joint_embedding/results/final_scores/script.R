cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)
library(rlang)

## VIASH START
par <- list(
  input = list.files("work/8a/76bbffc8075112714fba7fe14400f7", pattern = "*.h5ad$", full.names = TRUE),
  output = "output/pilot/joint_embedding/output.extract_scores.output.tsv",
  method_meta = NULL,
  metric_meta = list.files("src/joint_embedding/metrics", recursive = TRUE, pattern = "*.tsv$", full.names = TRUE),
  #solution_meta = "output/pilot/joint_embedding/meta_solution.collect_solution_metadata.output.tsv"
  dataset_meta = "results/meta_datasets.tsv"
)
## VIASH END

cat("Reading solution meta files\n")
dataset_meta <- readr::read_tsv(
  par$dataset_meta
)

cat("Reading metric meta files\n")
metric_defaults <-
  map_df(
    par$metric_meta,
    read_tsv,
    col_types = cols(
      metric_id = "c",
      metric_min = "c",
      metric_max = "c",
      metric_higherisbetter = "l"
    )
  ) %>%
  mutate(
    metric_min = as.numeric(metric_min),
    metric_max = as.numeric(metric_max)
  ) %>% 
  transmute(
    metric_id = metric_id, 
    missing_value = ifelse(metric_higherisbetter, metric_min, metric_max)
  )

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

  out <- as_tibble(ad$uns[c("dataset_id", "method_id", "metric_ids", "metric_values")])
  rm(ad)
  out
}) %>%
  rename(
    metric_id = metric_ids, 
    value = metric_values
  ) %>%
  mutate(
    dataset_orig_id = gsub("_JE$", "", dataset_id)
  )

expect_true(
  all(unique(scores$metric_id) %in% metric_defaults$metric_id),
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

cat("Creating default scores for missing entries based on metrics meta\n")
final_scores <- scores %>%
  full_join(method_meta, by = "method_id") %>%
  full_join(dataset_meta %>% select(dataset_orig_id = dataset_id, dataset_subtask = mod2_modality), by = "dataset_orig_id") %>%
  full_join(metric_defaults, by = "metric_id") %>%
  mutate(value_after_default = value %|% missing_value) %>%
  select(dataset_id, method_id, metric_id, dataset_subtask, value, value_after_default)

cat("Calculating geometric mean\n")
geomean <- final_scores %>%
  filter(metric_id %in% c("ari", "asw_batch", "asw_label", "cc_cons", "graph_conn", "nmi", "ti_cons_mean")) %>%
  group_by(dataset_id, method_id, dataset_subtask) %>%
  summarise_if(is.numeric, dynutils::calculate_geometric_mean) %>%
  ungroup() %>%
  mutate(metric_id = "geometric_mean")
final_scores <- bind_rows(final_scores %>% filter(metric_id != "geometric_mean"), geomean)

summary <-
  bind_rows(final_scores, final_scores %>% mutate(dataset_subtask = "Overall")) %>%
  group_by(method_id, metric_id, dataset_subtask) %>%
  summarise(
    mean = mean(value_after_default),
    var = var(value_after_default),
    .groups = "drop"
  )

# unique(scores$metric_id)

# summary %>%
#   filter(metric_id == "mean_spearman_per_cell") %>%
#   select(-var) %>%
#   spread(subtask, mean) %>%
#   arrange(Overall)

# summary %>%
#   filter(metric_id == "mse") %>%
#   select(-var) %>%
#   spread(subtask, mean) %>%
#   arrange(Overall)

jsontib <- summary %>%
  filter(metric_id == "geometric_mean") %>%
  select(-var) %>%
  spread(dataset_subtask, mean) %>%
  arrange(Overall)

cat("Writing output\n")
write_tsv(final_scores, par$output_scores)
write_tsv(summary, par$output_summary)
jsonlite::write_json(dynutils::tibble_as_list(jsontib), par$output_json)
