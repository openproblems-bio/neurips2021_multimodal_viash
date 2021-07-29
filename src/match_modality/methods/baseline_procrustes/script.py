import anndata
import scipy.spatial
import numpy as np

from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors

# VIASH START
par = {
    "input_mod1": "../../../../resources_test/match_modality/test_resource.mod1.h5ad",
    "input_mod2": "../../../../resources_test/match_modality/test_resource.mod2.h5ad",
    "output": "../../../../resources_test/match_modality/test_resource.prediction.h5ad",
    "n_svd": 100,
}
# VIASH END

# load dataset to be censored
ad_rna = anndata.read_h5ad(par["input_mod1"])
ad_mod2 = anndata.read_h5ad(par["input_mod2"])


rna, mod2 = ad_rna.X, ad_mod2.X

n_svd = min(par["n_svd"], rna.shape[1], mod2.shape[1])

rna_pca = TruncatedSVD(n_svd).fit_transform(rna)
mod2_pca = TruncatedSVD(n_svd).fit_transform(mod2)

rna_procrustes, mod2_procrustes, disparity = scipy.spatial.procrustes(rna_pca, mod2_pca)
print("Disparity value is: %0.3f" % disparity)

# 2. Perform nearest neighbors: find nearest neighbors for mod2 based on the rna
nn = NearestNeighbors(n_neighbors=1).fit(rna_procrustes)
distances, indices = nn.kneighbors(X=mod2_procrustes)

# 3. Helper: just a range -> so that each neighbor found with NN matches the right cell
indices_rna = list(range(rna.shape[0]))

pairing_matrix = np.zeros((rna.shape[0], mod2.shape[0]))
pairing_matrix[indices_rna, [x[0] for x in indices]] = 1

# Write out prediction
prediction = anndata.AnnData(
    X=pairing_matrix,
    uns={
        "dataset_id": ad_rna.uns["dataset_id"],
    }
)
prediction.write_h5ad(par["output"])
