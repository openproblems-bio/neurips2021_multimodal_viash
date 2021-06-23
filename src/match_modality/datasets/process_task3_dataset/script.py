# VIASH START
par = {
    "input": "resources/test/dyngen_bifurcating_antibody/dataset.h5ad",
    "output_censored": "output_censored.h5ad",
    "output_solution": "output_solution.h5ad"
}
# VIASH END

###############################################################################
###                            LOAD DEPENDENCIES                            ###
###############################################################################

import anndata
import random
import pandas as pd
import numpy as np

###############################################################################
###                             READ INPUT DATA                             ###
###############################################################################
# load dataset
adata = anndata.read_h5ad(par["input"])

# check which modalities are included
has_antibody = "protein" in adata.layers.keys()
has_atac = "chromatin" in adata.layers.keys()

# check if exactly 1 modality is included
# TODO: is this necessary? Can integration of all 3 modalities be required?
assert has_antibody != has_atac, "Strictly one of adata.layers[\"chromatin\"] and adata.layers[\"protein\"] must be " \
                                 "defined."

modality_types = { "modality1": "rna", "modality2": "chromatin" if has_atac else "protein"}

###############################################################################
###                          CREATE CENSOR OBJECT                           ###
###############################################################################

mod1 = adata.X
mod2 = adata.layers["protein"] if has_antibody else adata.layers["chromatin"]

# Generate the indices partition -> so that the shuffling can be saved
shuffle_cells = list(range(mod1.shape[0]))
random.shuffle(shuffle_cells)

# shuffle according to these indices
mod2 = mod2[shuffle_cells, :]

# create new anndata object
out_censor = anndata.AnnData(
    shape = adata.shape,
    layers = {
        "modality1": mod1,
        "modality2": mod2
    },
    uns = {
        "dataset_id": adata.uns["dataset_id"] + "_task3",
        "modality_types": modality_types
    },
    dtype = "float32"
)

out_censor.write_h5ad(filename = par["output_censored"], compression = "gzip")

###############################################################################
###                          CREATE SOLUTION OBJECT                         ###
###############################################################################

# create new anndata object
out_solution = anndata.AnnData(
    shape = adata.shape,
    obs = pd.DataFrame(
        index = adata.obs_names,
        data = { "shuffle_index": np.argsort(shuffle_cells) }
    ),
    uns = {
        "dataset_id": adata.uns["dataset_id"] + "_task3"
    }
)

out_solution.write_h5ad(filename = par["output_solution"], compression = "gzip")
