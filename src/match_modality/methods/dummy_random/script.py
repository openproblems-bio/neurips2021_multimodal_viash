import anndata as ad
import numpy as np
import scipy.sparse

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

# determine number of values in array
num_values = min(100, input_test_mod1.n_obs) * input_test_mod1.n_obs
indices = np.random.randint(input_test_mod1.n_obs**2, size=num_values)

mat_x = np.random.rand(num_values)
mat_i = indices % input_test_mod1.n_obs
mat_j = (indices / input_test_mod1.n_obs).astype(int)
pairing_matrix = scipy.sparse.csr_matrix(
    (mat_x, (mat_i, mat_j)),
    shape=(input_test_mod1.n_obs, input_test_mod2.n_obs)
)

# Write out prediction
prediction = ad.AnnData(
    X=pairing_matrix,
    uns={
        "method_id": "dummy_random",
        "dataset_id": input_test_mod1.uns["dataset_id"]
    }
)
prediction.write_h5ad(par["output"])
