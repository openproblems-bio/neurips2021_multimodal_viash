import anndata
import numpy as np

# VIASH START
par = {
    "input_mod1": "resources_test/match_modality/test_resource.mod1.h5ad",
    "input_mod2": "resources_test/match_modality/test_resource.mod2.h5ad",
    "output": "resources_test/match_modality/test_resource.prediction.h5ad",
}
# VIASH END

# load dataset to be censored
ad_mod1 = anndata.read_h5ad(par["input_mod1"])
ad_mod2 = anndata.read_h5ad(par["input_mod2"])

# Write out prediction
prediction = anndata.AnnData(
    X=np.ones((ad_mod1.n_obs, ad_mod2.n_obs)),
    uns={
        "method_id": "dummy_constant",
        "dataset_id": ad_mod1.uns["dataset_id"],
    }
)
prediction.write_h5ad(par["output"])
