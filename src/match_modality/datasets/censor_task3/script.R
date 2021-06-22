## VIASH START
par <- list(
  input = "../resources/test/dyngen_bifurcating_antibody/dataset.h5ad",
  # input = "resources/test/dyngen_bifurcating_atac/dataset.h5ad",
  output = "output.h5ad",
  shuffle_cells_output = "shuffle_cells.csv",
  shuffle_genes_output = "shuffle_genes.csv"
)
## VIASH END

# load libraries
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE)
library(readr)

# load h5ad file
adata <- anndata::read_h5ad(par$input)

# check for modalities
has_antibody <- "antibody" %in% names(adata$layers)
has_atac <- "atac" %in% names(adata$layers)

if (has_antibody == has_atac) {
  stop("Strictly one of adata.layers[\"atac\"] and adata.layers[\"antibody\"] must be defined.")
}

# hmm maybe not -> X always RNA?

rna <- adata$X
mod2 <- if (has_antibody) adata$layers[["antibody"]] else adata$layers[["atac"]]


# shuffle cells & genes
shuffle_cells <- sample.int(nrow(rna))
shuffle_genes <- sample.int(ncol(rna))

# TODO shuffle genen niet

rna <- rna[, shuffle_genes, drop = F]
mod2 <- mod2[shuffle_cells, shuffle_genes, drop = F]

write_csv(par$shuffle_cells_output, shuffle_cells)
write_csv(par$shuffle_genes_output, shuffle_genes)

# TODO: save this shuffling
# rename genes & remove cell names
rownames(rna) <- rownames(mod2) <- NULL
colnames(rna) <- colnames(mod2) <- paste0("gene_", seq_len(ncol(rna)))

# create anndata
censored <- anndata::AnnData(
  X = rna,
  shape = dim(counts),
  layers = list(
    modality2 = mod2
  ),
  uns = list(
    dataset_id = par$dataset_id,
    modality2 = if(has_antibody) "protein" else "atac"
  )
)

write_rds(out$model, par$output, compress = "gz")
