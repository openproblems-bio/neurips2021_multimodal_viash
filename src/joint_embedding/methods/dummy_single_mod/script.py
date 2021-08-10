import anndata
import numpy as np
import scanpy as sc

## VIASH START
par = {
    "input_mod1": "resources_test/joint_embedding/test_resource.mod1.h5ad",
    "input_mod2": "resources_test/joint_embedding/test_resource.mod2.h5ad",
    "output": "tmp/output_prediction.h5ad",
    "n_dims": 1,
    "modality": 1,
}
## VIASH END

assert 1 <= par["modality"] <= 2
print("Load and prepare data")

adata = anndata.read_h5ad(par["input_mod1" if par["modality"] == 1 else "input_mod2"])

sc.tl.pca(adata, n_comps=100)

print("Saving output")
adata_out = anndata.AnnData(
    X=adata.obsm["X_pca"],
    obs=adata.obs[["batch"]],
    uns={"dataset_id": adata.uns["dataset_id"], "method_id": "dummy_single_mod"},
)
adata_out.write_h5ad(par["output"], compression="gzip")
