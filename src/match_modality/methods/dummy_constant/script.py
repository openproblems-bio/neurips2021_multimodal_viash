import anndata as ad
import numpy as np

# VIASH START
par = {
    "input_train_mod1": "resources_test/match_modality/test_resource.train_mod1.h5ad",
    "input_train_mod2": "resources_test/match_modality/test_resource.train_mod2.h5ad",
    "input_train_sol": "resources_test/match_modality/test_resource.train_sol.h5ad",
    "input_test_mod1": "resources_test/match_modality/test_resource.test_mod1.h5ad",
    "input_test_mod2": "resources_test/match_modality/test_resource.test_mod2.h5ad",
    "output": "resources_test/match_modality/test_resource.prediction.h5ad",
}
# VIASH END

print("Load datasets")
input_test_mod1 = ad.read_h5ad(par["input_test_mod1"])
input_test_mod2 = ad.read_h5ad(par["input_test_mod2"])

print("Writing predictions")
prediction = ad.AnnData(
    X=np.ones((input_test_mod1.n_obs, input_test_mod2.n_obs)),
    uns={
        "method_id": "dummy_constant",
        "dataset_id": input_test_mod1.uns["dataset_id"],
    }
)
prediction.write_h5ad(par["output"])
