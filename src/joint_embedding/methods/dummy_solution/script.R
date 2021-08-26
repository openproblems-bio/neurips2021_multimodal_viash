cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = "resources_test/joint_embedding/test_resource.solution.h5ad"
)
## VIASH END

cat("Reading h5ad files\n")
ad_sol <- anndata::read_h5ad(par$input_solution)

ct <- ad_sol$obs$cell_type
pt1 <- ad_sol$obs$pseudotime_order_GEX
pt2 <- if (!is.null(ad_sol$obs$pseudotime_order_ADT)) ad_sol$obs$pseudotime_order_ADT else ad_sol$obs$pseudotime_order_ATAC

mat <- cbind(
  sapply(seq_len(100), function(i) {
    match(ct, sample(unique(ct)))
  }),
  sapply(seq_len(20), function(i) {
    ifelse(is.finite(pt1), pt1, runif(length(pt1)))
  }),
  sapply(seq_len(20), function(i) {
    ifelse(is.finite(pt2), pt2, runif(length(pt2)))
  })
)

cat("Performing DR\n")
dr <- dyndimred::dimred_knn_fr(mat, n_neighbors = round(mean(table(ct))))
rownames(dr) <- rownames(ad_sol)
colnames(dr) <- paste0("comp_", seq_len(ncol(dr)))

out <- anndata::AnnData(
  X = dr,
  uns = list(
    dataset_id = ad_sol$uns[["dataset_id"]],
    method_id = "dummy_solution"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
