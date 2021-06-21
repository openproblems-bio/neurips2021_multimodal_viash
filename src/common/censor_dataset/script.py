## VIASH START
par = {
  "input": "dataset.h5ad",
  "output": "dataset_censored.h5ad"
}

## VIASH END



# if a censored version of the dataset needs to be generated
# if (!is.null(par$output_censored)) {
#   counts <- dataset$X
#   counts_protein <- dataset$layers["counts_protein"]
#
#   # shuffle
#   shuffle_cells <- sample.int(nrow(counts_protein))
#   shuffle_genes <- sample.int(ncol(counts_protein))
#
#   counts <- counts[, shuffle_genes, drop = FALSE]
#   counts_protein <- counts_protein[shuffle_cells, shuffle_genes, drop = FALSE]
#
#   # rename genes, remove cell names
#   rownames(counts) <- rownames(counts_protein) <- NULL
#   colnames(counts) <- colnames(counts_protein) <- paste0("gene_", seq_len(ncol(counts)))
#
#   # create new dataset object
#   dataset_censored <- anndata::AnnData(
#     X = NULL,
#     shape = dim(counts),
#     layers = list(
#       modality1 = counts,
#       modality2 = counts_protein
#     ),
#     uns = list(
#       dataset_id = par$dataset_id,
#       modality1 = "mRNA",
#       modality2 = "protein"
#     )
#   )
#
#   dataset$uns$shuffle_cells <- shuffle_cells-1
#   dataset$uns$shuffle_genes <- shuffle_genes-1
#
#   dataset_censored$write_h5ad(par$output_censored, compression = "gzip")
# }
#
#
# dataset$write_h5ad(par$output, compression = "gzip")