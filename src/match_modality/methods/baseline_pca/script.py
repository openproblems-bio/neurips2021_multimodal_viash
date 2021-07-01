# VIASH START

par = {
    "input": "../../../../dataset_censored.h5ad",
    "output_rna": "output_rna.csv",
    "output_mod2": "output_mod2.csv"
}
# VIASH END

import anndata
import csv
import scipy.sparse
import numpy as np

from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

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
indices_paired = [(x, y[0]) for x, y in zip(indices_mod1, indices)]
pairing_matrix = np.zeros((rna.shape[0], mod2.shape[0]))

with open(par["output_rna"], 'w') as rna_file, open(par["output_mod2"], 'w') as mod2_file:
    rna_writer, mod2_writer = csv.writer(rna_file, delimiter=","), csv.writer(mod2_file, delimiter=",")
    rna_writer.writerows(rna_pca)
    mod2_writer.writerows(mod2_pca)
