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

pt_name <- if ("pseudotime_order_ATAC" %in% colnames(obs)) "pseudotime_order_ATAC" else "pseudotime_order_ADT"
obs_mat <- cbind(
  ct_mat,
  ifelse(is.na(obs[[pt_name]]), runif(nrow(obs)), obs[[pt_name]]),
  ifelse(is.na(obs$pseudotime_order_GEX), runif(nrow(obs)), obs$pseudotime_order_GEX)
)
dr <- uwot::umap(obs_mat, n_components = 10)

rownames(dr) <- rownames(ad_sol)
colnames(dr) <- paste0("comp_", seq_len(ncol(dr)))


# patchwork::wrap_plots(
#   map(colnames(obs), function(cn) {
#     vec <- obs[[cn]]
#     g <- 
#       qplot(dr[,1], dr[,2], colour = vec) + 
#       theme_bw() + 
#       labs(colour = cn)
#     if (is.numeric(vec)) {
#       g <- g + viridis::scale_colour_viridis()
#     }
#     g
#   })
# )

out <- anndata::AnnData(
  X = dr,
  uns = list(
    dataset_id = ad_sol$uns[["dataset_id"]],
    method_id = "dummy_solution"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
