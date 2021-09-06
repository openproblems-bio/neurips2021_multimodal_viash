cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = "output/datasets/joint_embedding/openproblems_bmmc_multiome_phase1/openproblems_bmmc_multiome_phase1.censor_dataset.output_solution.h5ad"
)
## VIASH END

cat("Reading h5ad files\n")
ad_sol <- anndata::read_h5ad(par$input_solution)
obs <- ad_sol$obs

ct_mat <- sapply(unique(obs$cell_type), function(ct) {
  (obs$cell_type == ct) + 0
})
obs_mat <- cbind(
  ct_mat,
  ifelse(is.na(obs$pseudotime_order_ATAC), runif(nrow(obs)), obs$pseudotime_order_ATAC),
  ifelse(is.na(obs$pseudotime_order_GEX), runif(nrow(obs)), obs$pseudotime_order_GEX)
)
dr <- prcomp(obs_mat)$x[,1:3]

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
