cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
out_path <- "output/pilot_inhouse/joint_embedding/output.final_scores.output_"
par <- list(
  input = list.files("/home/rcannood/workspace/viash_temp/neurips2021_work/85/c712592b7fd4839aef50eb4671f3f2/", pattern = "*.h5ad$", full.names = TRUE),
  method_meta = NULL,
  metric_meta = list.files("src/joint_embedding/metrics", recursive = TRUE, pattern = "*.tsv$", full.names = TRUE),
  dataset_meta = "output/datasets_2021-11-08/phase2_private/meta.tsv",
  output_scores = paste0(out_path, "scores.tsv"),
  output_summary = paste0(out_path, "summary.tsv"),
  output_json = paste0(out_path, "json.json")
)
## VIASH END

json_metrics <- c("asw_batch", "asw_label", "cc_cons", "graph_conn", "nmi", "ti_cons_batch_mean", "arithmetic_mean")

cat("Reading solution meta files\n")
dataset_meta <- 
  readr::read_tsv(par$dataset_meta) %>% 
  transmute(
    dataset_orig_id = dataset_id,
    dataset_id = paste0(dataset_id, "_JE"),
    dataset_subtask = mod2_modality
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
default_scores <- crossing(
  method_id = method_meta$method_id, 
  dataset_meta %>% select(dataset_id, dataset_orig_id),
  metric_defaults %>% mutate(value = NA_real_, value_after_default = missing_value) %>% select(-missing_value)
)

scores1 <- bind_rows(
  scores %>% mutate(value_after_default = value),
  anti_join(default_scores, scores, by = c("method_id", "dataset_id", "metric_id"))
) %>% 
  left_join(dataset_meta, by = c("dataset_id", "dataset_orig_id")) %>%
  filter(metric_id != "latent_mixing", !grepl("rfoob_", metric_id))

cat("Calculating geometric mean\n")
arimean <- scores1 %>%
  filter(metric_id %in% c("asw_batch", "asw_label", "cc_cons", "graph_conn", "nmi", "ti_cons_batch_mean")) %>%
  group_by(dataset_id, dataset_orig_id, method_id) %>%
  summarise_if(is.numeric, function(x) dynutils::calculate_arithmetic_mean(x)) %>%
  ungroup() %>%
  mutate(metric_id = "arithmetic_mean")
final_scores <- bind_rows(scores1 %>% filter(metric_id != "arithmetic_mean"), arimean)

summary <-
  final_scores %>%
  # bind_rows(final_scores, final_scores %>% mutate(dataset_subtask = "Overall")) %>%
  group_by(method_id, metric_id, dataset_subtask) %>%
  summarise(
    mean = mean(value_after_default),
    var = var(value_after_default),
    .groups = "drop"
  )

# unique(scores$metric_id)

json_out <- summary %>%
  filter(metric_id %in% json_metrics) %>%
  mutate(comb_id = ifelse(metric_id == "arithmetic_mean", metric_id, paste0(metric_id, "_", dataset_subtask))) %>%
  select(-var, -metric_id, -dataset_subtask) %>%
  spread(comb_id, mean) %>%
  arrange(asw_label_ADT)

cat("Writing output\n")
final_scores <- final_scores %>% map(as.vector) %>% as_tibble
summary <- summary %>% map(as.vector) %>% as_tibble
readr::write_tsv(final_scores, par$output_scores)
readr::write_tsv(summary, par$output_summary)
jsonlite::write_json(json_out, par$output_json)
