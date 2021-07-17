cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input = list.files("/tmp/neurips2021_work", full.names = TRUE, recursive = TRUE, pattern = "*calculate_task1_metrics.output.h5ad"),
  output = "output/task1_scores.tsv",
  summary = "output/task1_summary.tsv"
)
par$input <- par$input[!duplicated(basename(par$input))]
inp <- par$input[[1]]
## VIASH END


cat("Reading input h5ad files")
scores <- map_df(par$input, function(inp) {
  cat("Reading '", inp, "'\n", sep = "")
  ad <- read_h5ad(inp)

  for (uns_name in c("dataset_id", "method_id", "metric_ids", "metric_values")) {
    expect_true(
      uns_name %in% names(ad$uns),
      info = paste0("File ", inp, " must contain `uns['", uns_name, "']`")
    )
  }

  as_tibble(ad$uns[c("dataset_id", "method_id", "metric_ids", "metric_values")])
})

summary <- 
  scores %>% 
  group_by(method_id, metric_ids) %>% 
  summarise(metric_values = mean(metric_values), .groups = "drop")

score_spread <- scores %>% spread(metric_ids, metric_values)
summary_spread <- summary %>% spread(metric_ids, metric_values)

write_tsv(score_spread, par$output)
write_tsv(summary_spread, par$summary)
