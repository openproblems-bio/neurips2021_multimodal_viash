# Dependencies:
# pip: scikit-learn, anndata, scanpy
#
# Python starter kit for the NeurIPS 2021 Single-Cell Competition.
# Parts with `TODO` are supposed to be changed by you.
#
# More documentation:
#
# https://viash.io/docs/creating_components/python/

import logging
import anndata as ad
import numpy as np
import scipy

from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
logging.basicConfig(level=logging.INFO)

## VIASH START
# Anything within this block will be removed by `viash` and will be
# replaced with the parameters as specified in your config.vsh.yaml.

dataset_path = 'sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.'
# dataset_path = 'output/datasets/match_modality/openproblems_bmmc_multiome_phase1_rna/openproblems_bmmc_multiome_phase1_rna.censor_dataset.output_'

par = {
    'input_train_mod1': f'{dataset_path}train_mod1.h5ad',
    'input_train_mod2': f'{dataset_path}train_mod2.h5ad',
    'input_train_sol': f'{dataset_path}train_sol.h5ad',
    'input_test_mod1': f'{dataset_path}test_mod1.h5ad',
    'input_test_mod2': f'{dataset_path}test_mod2.h5ad',
    'distance_method': 'minkowski',
    'output': 'output.h5ad',
    'n_pcs': 4,
    'n_neighbors': 5,
}
## VIASH END

# TODO: change this to the name of your method
method_id = "python_starter_kit"

logging.info('Reading `h5ad` files...')
input_train_mod1 = ad.read_h5ad(par['input_train_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])
input_train_sol = ad.read_h5ad(par['input_train_sol'])
input_test_mod1 = ad.read_h5ad(par['input_test_mod1'])
input_test_mod2 = ad.read_h5ad(par['input_test_mod2'])

# unscramble training cells
ord = np.argsort(input_train_sol.uns['pairing_ix'])
input_train_mod2 = input_train_mod2[ord, :]


# TODO: implement own method

# This starter kit is split up into several steps.
# * compute dimensionality reduction on [train_mod1, test_mod1] data
# * train linear model to predict the train_mod2 data from the dr_mod1 values
# * predict test_mod2 matrix from model and test_mod1
# * calculate k nearest neighbors between test_mod2 and predicted test_mod2
# * transform k nearest neighbors into a pairing matrix

# TODO: implement own method
input_train = ad.concat(
    {"train": input_train_mod1, "test": input_test_mod1},
    axis=0,
    join="outer",
    label="group",
    fill_value=0,
    index_unique="-"
)

# Do PCA on the input data
logging.info('Performing dimensionality reduction on modality 1 values...')
embedder_mod1 = TruncatedSVD(n_components=50)
pca_mod1 = embedder_mod1.fit_transform(input_train.X)

logging.info('Performing dimensionality reduction on modality 2 values...')
embedder_mod2 = TruncatedSVD(n_components=50)
# input_train_mod2$X[order(match_train)
pca_train_mod2 = embedder_mod2.fit_transform(input_train_mod2.X)

# split dimred back up
pca_train_mod1 = pca_mod1[input_train.obs['group'] == 'train']
pca_test_mod1 = pca_mod1[input_train.obs['group'] == 'test']

assert len(pca_train_mod1) + len(pca_test_mod1) == len(pca_mod1)

logging.info('Running Linear regression...')
reg = LinearRegression()

# Train the model on the PCA reduced modality 1 and 2 data
reg.fit(pca_train_mod1, pca_train_mod2)
pca_test_mod2_pred = reg.predict(pca_test_mod1)

# Project the predictions back to the modality 2 feature space
input_test_mod2_pred = pca_test_mod2_pred @ embedder_mod2.components_


logging.info('calculate k nearest neighbors between test_mod2 and predicted test_mod2')
# To get the matching matrix, for each point in mod1_test, we take the 1000 nearest neighbors of that
# point in the transformed mod2_test dataset
n_neighbors = min(par["n_neighbors"], input_test_mod1.n_obs)
nn = NearestNeighbors(n_neighbors=n_neighbors).fit(input_test_mod2_pred)
distances, indices = nn.kneighbors(X=input_test_mod2.X)

logging.info('transform k nearest neighbors into a pairing matrix')
# Translate the neighborhood assignments to a pairing matrix that is (n_obs, n_obs)
# NOTE: `pairing_matrix` must have NO MORE than 1000*n_obs non-zero entries
ind_i = np.tile(np.arange(input_test_mod1.n_obs), (n_neighbors, 1)).T.flatten()
ind_j = indices.flatten()
ind_dist = distances.flatten()
ind_x = 2 * max(ind_dist) - ind_dist
pairing_matrix = scipy.sparse.csr_matrix(
    (ind_x, (ind_i, ind_j)),
    shape=(input_test_mod1.n_obs, input_test_mod2.n_obs)
)

# Normalize values to sum to 1
pairing_matrix = normalize(pairing_matrix, norm="l1", axis=1)


logging.info('write prediction output')
out = ad.AnnData(
    X=pairing_matrix,
    uns={
        "dataset_id": input_train_mod1.uns["dataset_id"],
        "method_id": method_id
    }
)
out.write_h5ad(par['output'], compression="gzip")
