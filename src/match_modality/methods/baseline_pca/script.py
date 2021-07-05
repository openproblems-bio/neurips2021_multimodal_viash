import anndata
import scipy.sparse
import numpy as np

from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

# VIASH START
par = {
    "input": "../../../../dataset_censored.h5ad",
    "prediction": "prediction.h5ad"
}
# VIASH END

# load dataset to be censored
adata = anndata.read_h5ad(par["input"])

rna = adata.X
mod2 = adata.layers["modality2"]

comb = scipy.sparse.vstack([rna, mod2]).transpose().todense()

pca = PCA(n_components=10)
pca.fit(comb)

# find kNN in the PCA

rna_pca = pca.components_[:, :rna.shape[0]].transpose()
mod2_pca = pca.components_[:, rna.shape[0]:].transpose()

nn = NearestNeighbors(n_neighbors=1).fit(rna_pca)
distances, indices = nn.kneighbors(X=mod2_pca)
indices_mod1 = [i for i in range(rna.shape[0])]

# indices_paired = [(x, np.array([y[0], z[0]])) for x, y, z in zip(indices_mod1, indices, distances)]
# indices_paired2 = [[x, y[0]] for x, y in zip(indices_mod1, indices)]
# adata.obsp = indices_paired

pairing_matrix = np.zeros((rna.shape[0], mod2.shape[0]))
# pairing_matrix[indices_paired2] = 1

# TODO coo matrix?
pairing_matrix[indices_mod1, [x[0] for x in indices]] = 1

adata.obsp["result"] = pairing_matrix
adata.write_h5ad(par["prediction"])

