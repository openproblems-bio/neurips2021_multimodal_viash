# LOAD DEPENDENCIES

import anndata
import random
import pandas as pd
import numpy as np
import scipy.sparse

# VIASH START
par = {
    "input_rna": "pbmc_1k_protein_v3.output_rna.h5ad",
    "input_mod2": "pbmc_1k_protein_v3.output_mod2.h5ad",
    "output_mod1": "pbmc_1k_protein_v3.censored_rna.h5ad",
    "output_mod2": "pbmc_1k_protein_v3.censored_mod2.h5ad",
    "output_solution": "pbmc_1k_protein_v3.solution.h5ad",
}
# VIASH END

# READ INPUT DATA
rna = anndata.read_h5ad(par["input_rna"])
mod2 = anndata.read_h5ad(par["input_mod2"])

# TODO: change on how to know which modality is there
# modality_types = {"modality1": "rna", "modality2": "chromatin" if has_atac else "protein"}

#  CREATE CENSOR OBJECT
pairings = scipy.sparse.spdiags(np.full(rna.shape[0], 1), diags=0, m=rna.shape[0], n=mod2.shape[0], format="csr")

# Generate the indices partition
shuffle_cells = list(range(rna.shape[0]))
random.shuffle(shuffle_cells)

# shuffle according to these indices
mod2 = mod2[shuffle_cells, :]
pairings = pairings[:, shuffle_cells]

censor_rna = anndata.AnnData(
    X=rna.X,
    uns={
        "dataset_id": rna.uns["dataset_id"] + "_task3",
        "modality": "rna",
        "paired": True,
    },
    dtype="float32",
)

censor_mod2 = anndata.AnnData(
    X=mod2.X,
    uns={
        "dataset_id": mod2.uns["dataset_id"] + "_task3",
        "modality": mod2.uns["modality"],
        "paired": True,
    },
    dtype="float32",
)

censor_rna.write_h5ad(filename=par["output_mod1"], compression="gzip")
censor_mod2.write_h5ad(filename=par["output_mod2"], compression="gzip")

# CREATE SOLUTION PAIRING
out_solution = anndata.AnnData(
    shape=rna.shape,
    uns={
        "dataset_id": rna.uns["dataset_id"] + "_task3",
        "pairings": pairings,
    },
    dtype="float32"
)

out_solution.write_h5ad(filename=par["output_solution"], compression="gzip")
