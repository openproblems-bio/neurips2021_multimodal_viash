# VIASH START
par = {
    "input": "../../../resources/test/dyngen_bifurcating_antibody/dataset.h5ad",
    "output": "output.h5ad",
    "shuffle_cells_output": "shuffle_cells.csv",
    "shuffle_genes_output": "shuffle_genes.csv"
}
# VIASH END

import anndata
import random

# load dataset to be censored
adata = anndata.read_h5ad(par["input"])

# check which modalities are included
has_antibody = "antibody" in adata.layers.keys()
has_atac = "atac" in adata.layers.keys()

# check if exactly 1 modality is included
# TODO: is this necessary? Can integration of all 3 modalities be required?
assert has_antibody != has_atac, "Strictly one of adata.layers[\"atac\"] and adata.layers[\"antibody\"] must be " \
                                 "defined."

rna = adata.X
mod2 = adata.layers["antibody"] if has_antibody else adata.layers["atac"]

# Generate the indices partition -> so that the shuffling can be saved
shuffle_cells = random.sample([x for x in range(rna.shape[0])], k=rna.shape[0])
shuffle_genes = random.sample([x for x in range(rna.shape[1])], k=rna.shape[1])

# shuffle according to these indices
rna = rna[:, shuffle_genes]
mod2 = mod2[:, shuffle_genes]
mod2 = mod2[shuffle_cells, :]

# TODO: write shuffling to a csv

# genes & cell names -> need to change gene names & save in anndata, as cell names are not present in the extracted data

# write_csv(par$shuffle_cells_output, shuffle_cells)
# write_csv(par$shuffle_genes_output, shuffle_genes)

censored = anndata.AnnData(
    X=rna,
    shape=rna.shape,
    layers={
        "modality2": mod2
    },
    uns={
        "dataset_id":par["dataset_id"],
        "modality2": "protein" if has_antibody else "atac"
    }
)

# # TODO: save this shuffling
# # rename genes & remove cell names
# rownames(rna) <- rownames(mod2) <- NULL
# colnames(rna) <- colnames(mod2) <- paste0("gene_", seq_len(ncol(rna)))
#
# # create anndata
# censored <- anndata::AnnData(
#   X = rna,
#   shape = dim(counts),
#   layers = list(
#     modality2 = mod2
#   ),
#   uns = list(
#     dataset_id = par$dataset_id,
#     modality2 = if(has_antibody) "protein" else "atac"
#   )
# )
#
# write_rds(out$model, par$output, compress = "gz")
