import anndata
from scipy.sparse import csr_matrix
import numpy as np

# VIASH START
par = {
    "input_train_mod1": "resources_test/predict_modality/test_resource.train_mod1.h5ad",
    "input_test_mod1": "resources_test/predict_modality/test_resource.test_mod1.h5ad",
    "input_train_mod2": "resources_test/predict_modality/test_resource.train_mod2.h5ad",
    "output": "output.h5ad",
}
# VIASH END

# load dataset to be censored
ad_mod1_train = anndata.read_h5ad(par["input_train_mod1"])
ad_mod1_test = anndata.read_h5ad(par["input_test_mod1"])
ad_mod2 = anndata.read_h5ad(par["input_train_mod2"])

# Testing with sparse prediction matrix
prediction = csr_matrix((ad_mod1_test.n_obs, ad_mod2.n_vars), dtype = np.float32)

# Write out prediction
out = anndata.AnnData(
    X=prediction,
    uns={
        "dataset_id": ad_mod2.uns["dataset_id"],
        "method_id": "dummy_zeros",
    }
)
out.write_h5ad(par["output"])
