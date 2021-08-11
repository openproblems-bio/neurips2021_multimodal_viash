import anndata
from scipy.sparse import csr_matrix
import numpy as np

# VIASH START
par = {
    "input_mod1_train": "resources_test/predict_modality/test_resource.mod1.train.h5ad",
    "input_mod1_test": "resources_test/predict_modality/test_resource.mod1.test.h5ad",
    "input_mod2": "resources_test/predict_modality/test_resource.mod2.h5ad",
    "output": "test_resource.prediction.h5ad",
}
# VIASH END

# load dataset to be censored
ad_mod1_train = anndata.read_h5ad(par["input_mod1_train"])
ad_mod1_test = anndata.read_h5ad(par["input_mod1_test"])
ad_mod2 = anndata.read_h5ad(par["input_mod2"])
ad_test = ad_mod1[ad_mod1.obs["group"] == "test"]

# Testing with sparse prediction matrix
prediction = csr_matrix((ad_test.n_obs, ad_mod2.n_vars), dtype = np.float32)

# Write out prediction
out = anndata.AnnData(
    X=prediction,
    uns={
        "dataset_id": ad_mod1.uns["dataset_id"],
        "method_id": "dummy_zeros",
    }
)
out.write_h5ad(par["output"])
