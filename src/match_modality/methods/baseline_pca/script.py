import anndata
import scipy.sparse
import numpy as np

from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors

# VIASH START
par = {
    "input_mod1": "../../../../resources_test/common/pbmc_1k_protein_v3.censored_rna.h5ad",
    "input_mod2": "../../../../resources_test/common/pbmc_1k_protein_v3.censored_mod2.h5ad",
    "output": "../../../../resources_test/common/pbmc_1k_protein_v3.prediction.h5ad",
}
# VIASH END

# load dataset to be censored
ad_rna = anndata.read_h5ad(par["input_mod1"])
ad_mod2 = anndata.read_h5ad(par["input_mod2"])

rna, mod2 = ad_rna.X, ad_mod2.X

rna_y = rna.shape[1]
mod2_y = mod2.shape[1]
max_y = max(rna_y, mod2_y)

# pad datasets with zeroes, so that they contain the same amount of features
# Necessary for performing the truncated SVD
rna.resize(rna.shape[0], max_y)
mod2.resize(mod2.shape[0], max_y)

# Shape of comb should be: (n_cells, n_features)
comb = scipy.sparse.vstack([rna, mod2])
tsvd = TruncatedSVD(n_components=10)
dimred = tsvd.fit_transform(comb)

# find the nearest neighbors accross modalities in the PCA
# 1. split the dimred into rna & mod2 cells
rna_pca = dimred[:rna.shape[0], :]
mod2_pca = dimred[rna.shape[0]:, :]

# 2. Perform nearest neighbors: find nearest neighbors for mod2 based on the rna
nn = NearestNeighbors(n_neighbors=1).fit(rna_pca)
distances, indices = nn.kneighbors(X=mod2_pca)

# 3. Helper: just a range -> so that each neighbor found with NN matches the right cell
indices_rna = [i for i in range(rna.shape[0])]

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
