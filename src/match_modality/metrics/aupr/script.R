cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("pracma", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = "resources_test/match_modality/test_resource.test_sol.h5ad",
  input_prediction = "resources_test/match_modality/test_resource.prediction.h5ad",
  output = "resources_test/match_modality/test_resource.scores.h5ad"
)
## VIASH END

cat("Read solution h5ad\n")
ad_sol <- anndata::read_h5ad(par$input_solution)

cat("Read prediction h5ad\n")
expect_true(
  grepl("\\.h5ad$", par$input_prediction),
  info = "Prediction file should be an h5ad file"
)
ad_pred <-
  tryCatch({
    anndata::read_h5ad(par$input_prediction)
  }, error = function(e) {
    stop(paste0("Can't open prediction h5ad file. Detailed error message:\n", e$message))
  })
expect_true(
  ad_sol$uns$dataset_id == ad_pred$uns$dataset_id
)
X_sol <- ad_sol$X
X_pred <- as(ad_pred$X, "CsparseMatrix")
dimnames(X_sol) <- dimnames(X_pred) <- list(NULL, NULL)

cell_type <- ad_sol$obs$cell_type

cat("Data wrangling\n")
sol_summ <- summary(X_sol) %>%
  as_tibble() %>%
  rename(gold = x) %>%
  filter(gold != 0)
pred_summ <- summary(X_pred) %>%
  left_join(sol_summ, by = c("i", "j")) %>%
  as_tibble() %>%
  mutate(
    gold = ifelse(is.na(gold), 0, gold),
    label_match = (cell_type[i] == cell_type[ad_sol$uns$pairing_ix+1][j])+0
  )

expect_true(
  nrow(pred_summ) <= 1000 * nrow(sol_summ),
  info = "Number of non-zero values for the prediction should be less or equal to 1000 times the number of cells in the dataset."
)

# helper function
calculate_au <- function(values, are_true, num_positive_interactions, num_possible_interactions, extend_by = 10000) {
  ord <- order(rank(values, ties.method = "random"), decreasing = TRUE)
  values <- values[ord]
  are_true <- are_true[ord]

  # calculate base statistics
  num_selected <- seq_along(are_true)
  tp <- cumsum(are_true)
  fp <- num_selected - tp
  length_ranking <- length(tp)
  num_negative_interactions <- num_possible_interactions - num_positive_interactions

  # extend base statistics, if necessary
  if (extend_by > 0 && length_ranking != num_possible_interactions) {
    diff.predictions <- num_possible_interactions - length_ranking
    diff.trues <- num_positive_interactions - tail(tp, 1)
    diff.negs <- num_negative_interactions - tail(fp, 1)

    multiplier <- seq_len(extend_by) / extend_by

    extra_num_selected <- multiplier * diff.predictions + tail(num_selected, 1)
    extra_tp <- multiplier * diff.trues + tail(tp, 1)
    extra_fp <- multiplier * diff.negs + tail(fp, 1)

    num_selected <- c(num_selected, extra_num_selected)
    are_true <- c(are_true, rep(NA, extend_by))
    tp <- c(tp, extra_tp)
    fp <- c(fp, extra_fp)
  }

  # calculate extended statistics
  metrics <- tibble(
    num_selected = c(0, num_selected),
    are_true = c(NA, are_true),
    tp = c(0, tp),
    fp = c(0, fp),
    fn = num_positive_interactions - tp,
    tn = num_negative_interactions - fp,
    acc = (tp + tn) / (num_positive_interactions + num_negative_interactions),
    tpr = tp / num_positive_interactions,
    spec = tn / num_negative_interactions,
    prec = ifelse(num_selected == 0, 1, tp / (tp + fp)),
    npv = tn / (tn + fn),
    f1 = 2 * tp / (2 * tp + fp + fn),
    mcc = ifelse(num_selected == 0, 0, (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))),
    informedness = tpr + spec - 1,
    markedness = prec + npv - 1
  )

  # calculate area under the curves
  area_under <- tibble(
    auroc = pracma::trapz(1 - metrics$spec, metrics$tpr),
    aupr = abs(pracma::trapz(metrics$tpr, metrics$prec))
  )

  list(metrics = metrics, area_under = area_under)
}


cat("Calculate area under the curve\n")
au_out <- calculate_au(
  values = pred_summ$x,
  are_true = pred_summ$gold,
  num_positive_interactions = nrow(sol_summ),
  num_possible_interactions = (nrow(ad_sol) * 1.0) * nrow(ad_sol)
)

label_numbers <- as.numeric(table(cell_type))
au_match_out <- calculate_au(
  values = pred_summ$gold,
  are_true = pred_summ$label_match,
  num_positive_interactions = sum((label_numbers - 1) * label_numbers),
  num_possible_interactions = (nrow(ad_sol) * 1.0) * nrow(ad_sol)
)

# GENIE3bis::plot_curves(au_out)
# GENIE3bis::plot_curves(au_match_out)

colnames(au_out$area_under) <- paste0("pairing_", colnames(au_out$area_under))
colnames(au_match_out$area_under) <- paste0("celltype_", colnames(au_match_out$area_under))

cat("Create output object\n")
out_values <- c(as.list(au_out$area_under), as.list(au_match_out$area_under))

out <- anndata::AnnData(
  X = NULL,
  shape = dim(ad_sol),
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = names(out_values),
    metric_values = as.numeric(out_values)
  )
)

# should we also save the metrics object?
# this would allow for plotting the auroc and aupr curves afterwards.

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")