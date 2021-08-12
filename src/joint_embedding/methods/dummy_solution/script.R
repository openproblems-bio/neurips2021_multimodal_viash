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

mat <- sapply(seq_len(100), function(i) {
  match(ct, sample(unique(ct)))
})

cat("Performing DR\n")
dr <- lmds::lmds(mat, ndim = 10, distance_method = "pearson")

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
