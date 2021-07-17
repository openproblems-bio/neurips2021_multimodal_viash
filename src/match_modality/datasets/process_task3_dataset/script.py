print("Loading dependencies")
import anndata
import random
import pandas as pd
import numpy as np
import scipy.sparse

# VIASH START
par = {
    "input_mod1": "resources_test/common/test_resource.output_rna.h5ad",
    "input_mod2": "resources_test/common/test_resource.output_mod2.h5ad",
    "output_mod1": "resources_test/task3/test_resource.mod1.h5ad",
    "output_mod2": "resources_test/task3/test_resource.mod2.h5ad",
    "output_solution": "resources_test/task3/test_resource.solution.h5ad",
}
# VIASH END

print("Reading input data")
mod1 = anndata.read_h5ad(par["input_mod1"])
mod2 = anndata.read_h5ad(par["input_mod2"])

print("Shuffling rows of mod2")
pairings = scipy.sparse.spdiags(np.full(mod1.shape[0], 1), diags=0, m=mod1.shape[0], n=mod2.shape[0], format="csr")

# Generate the indices partition
shuffle_cells = list(range(mod1.shape[0]))
random.shuffle(shuffle_cells)

print("Creating output objects")
desired_var1_cols = [x for x in ["gene_ids", "feature_types"] if x in mod1.var.columns]
out_mod1 = anndata.AnnData(
    X=mod1.X,
    var=mod1.var[desired_var1_cols],
    uns={
        "dataset_id": mod1.uns["dataset_id"] + "_task3"
    },
    dtype="float32",
)
out_mod1.X.sort_indices()

desired_var2_cols = [x for x in ["gene_ids", "feature_types"] if x in mod2.var.columns]
out_mod2 = anndata.AnnData(
    X=mod2.X[shuffle_cells, :],
    var=mod2.var[desired_var2_cols],
    uns={
        "dataset_id": mod2.uns["dataset_id"] + "_task3"
    },
    dtype="float32",
)
out_mod2.X.sort_indices()

out_solution = anndata.AnnData(
    X=pairings[:, shuffle_cells],
    uns={
        "dataset_id": mod1.uns["dataset_id"] + "_task3"
    },
    dtype="float32"
)

print("Writing output objects to file")
out_mod1.write_h5ad(filename=par["output_mod1"], compression="gzip")
out_mod2.write_h5ad(filename=par["output_mod2"], compression="gzip")
out_solution.write_h5ad(filename=par["output_solution"], compression="gzip")
