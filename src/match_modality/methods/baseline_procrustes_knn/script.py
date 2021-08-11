import anndata as ad
import scipy.spatial
import scipy.sparse
import numpy as np

from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors

# VIASH START
par = {
    "input_train_mod1": "resources_test/match_modality/test_resource.train_mod1.h5ad",
    "input_train_mod2": "resources_test/match_modality/test_resource.train_mod2.h5ad",
    "input_train_sol": "resources_test/match_modality/test_resource.train_sol.h5ad",
    "input_test_mod1": "resources_test/match_modality/test_resource.test_mod1.h5ad",
    "input_test_mod2": "resources_test/match_modality/test_resource.test_mod2.h5ad",
    "output": "resources_test/match_modality/test_resource.prediction.h5ad",
    "n_svd": 100,
}
# VIASH END

print("Load datasets")
input_train_mod1 = ad.read_h5ad(par["input_train_mod1"])
input_train_mod2 = ad.read_h5ad(par["input_train_mod2"])
# input_train_sol = ad.read_h5ad(par["input_train_sol"])
input_test_mod1 = ad.read_h5ad(par["input_test_mod1"])
input_test_mod2 = ad.read_h5ad(par["input_test_mod2"])

# concatenate train and test data
mod1 = ad.concat(
    { "train": input_train_mod1, "test": input_test_mod1 }, 
    index_unique="-",
    label="group"
)
mod2 = ad.concat(
    { "train": input_train_mod2, "test": input_test_mod2 }, 
    index_unique="-",
    label="group"
)
# create helper views
mod1te = mod1[mod1.obs["group"] == "test", :]
mod2te = mod2[mod2.obs["group"] == "test", :]

print("Perform PCA")
n_svd = min(par["n_svd"], mod1.n_obs, mod2.n_obs, mod1.n_vars, mod1.n_vars)

mod1.obsm["X_pca"] = TruncatedSVD(n_svd).fit_transform(mod1.X)
mod2.obsm["X_pca"] = TruncatedSVD(n_svd).fit_transform(mod2.X)

print("Run procrustes")
mod1.obsm["X_pro"], mod2.obsm["X_pro"], disparity = scipy.spatial.procrustes(
    mod1.obsm["X_pca"], 
    mod2.obsm["X_pca"]
)
print("> Disparity value is: %0.3f" % disparity)

print("Perform nearest neighbors")
n_neighbors = min(100, mod1te.n_obs)
nn = NearestNeighbors(n_neighbors=n_neighbors).fit(mod1te.obsm["X_pro"])
distances, indices = nn.kneighbors(X=mod2te.obsm["X_pro"])

print("Create pairing matrix")
ind_i = np.tile(np.arange(mod1te.n_obs), (n_neighbors, 1)).T.flatten()
ind_j = indices.flatten()
ind_dist = distances.flatten()
ind_x = 2 * max(ind_dist) - ind_dist
pairing_matrix = scipy.sparse.csr_matrix((ind_x, (ind_i, ind_j)))

print("Write prediction output")
prediction = ad.AnnData(
    X=pairing_matrix,
    uns={
        "dataset_id": input_train_mod1.uns["dataset_id"],
        "method_id": "baseline_procrustes_knn"
    }
)
prediction.write_h5ad(par["output"])
