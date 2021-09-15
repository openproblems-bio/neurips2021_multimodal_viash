cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
out_path <- "output/pilot_inhouse/match_modality/output.final_scores.output_"
par <- list(
  input = list.files("/tmp/neurips2021_work/d0/a319268c0a70d0890301ad3dd8628d/", pattern = "*.h5ad$", full.names = TRUE),
  method_meta = NULL,
  metric_meta = list.files("src/match_modality/metrics", recursive = TRUE, pattern = "*.tsv$", full.names = TRUE),
  #solution_meta = "output/pilot/match_modality/meta_solution.collect_solution_metadata.output.tsv"
  dataset_meta = "results/meta_datasets.tsv",
  output_scores = paste0(out_path, "scores.tsv"),
  output_summary = paste0(out_path, "summary.tsv"),
  output_json = paste0(out_path, "json.json")
)
## VIASH END

json_metric <- "match_probability"

cat("Reading solution meta files\n")
dm <- readr::read_tsv(par$dataset_meta)
dataset_meta <-
  bind_rows(dm, dm %>% rename(mod1_modality = mod2_modality, mod2_modality = mod1_modality)) %>% 
  transmute(
    dataset_orig_id = dataset_id,
    dataset_subtask = paste0(mod1_modality, "2", mod2_modality),
    dataset_id = paste0(dataset_orig_id, "_MM_", tolower(dataset_subtask))
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
    dataset_orig_id = gsub("_MM.*$", "", dataset_id)
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
default_scores <- crossing(
  method_id = method_meta$method_id, 
  dataset_meta %>% select(dataset_id, dataset_orig_id),
  metric_defaults %>% mutate(value = NA_real_, value_after_default = missing_value) %>% select(-missing_value)
)

final_scores <- bind_rows(
  scores %>% mutate(value_after_default = value),
  anti_join(default_scores, scores, by = c("method_id", "dataset_id", "metric_id"))
) %>% 
  left_join(dataset_meta, by = c("dataset_id", "dataset_orig_id"))

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
#   filter(metric_id == "match_probability") %>%
#   select(-var) %>%
#   spread(dataset_subtask, mean) %>%
#   arrange(Overall)

json_out <- summary %>%
  filter(metric_id == json_metric) %>%
  select(-var, -metric_id) %>%
  spread(dataset_subtask, mean) %>%
  arrange(Overall)

cat("Writing output\n")
final_scores <- final_scores %>% map(as.vector) %>% as_tibble
summary <- summary %>% map(as.vector) %>% as_tibble
write_tsv(final_scores, par$output_scores)
write_tsv(summary, par$output_summary)
jsonlite::write_json(json_out, par$output_json)