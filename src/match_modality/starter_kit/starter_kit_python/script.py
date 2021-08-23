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
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neighbors import NearestNeighbors

logging.basicConfig(level=logging.INFO)

## VIASH START
# Anything within this block will be removed by `viash` and will be
# replaced with the parameters as specified in your config.vsh.yaml.
par = {
    'input_train_mod1': 'resources_test/match_modality/test_resource.train_mod1.h5ad',
    'input_train_mod2': 'resources_test/match_modality/test_resource.train_mod2.h5ad',
    'input_train_sol': 'resources_test/match_modality/test_resource.train_sol.h5ad',
    'input_test_mod1': 'resources_test/match_modality/test_resource.test_mod1.h5ad',
    'input_test_mod2': 'resources_test/match_modality/test_resource.test_mod2.h5ad',
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


# TODO: implement own method

# This starter kit is split up into several steps.
# * compute dimensionality reduction on [train_mod1, test_mod1] data
# * train regression model to predict the train_mod2 data from the dr_mod1 values
# * predict test_mod2 matrix from model and test_mod1
# * calculate k nearest neighbors between test_mod2 and predicted test_mod2
# * transform k nearest neighbors into a pairing matrix

logging.info('Compute dimensionality reduction on [train_mod1, test_mod1] data')
input_train = ad.concat(
    {"train": input_train_mod1, "test": input_test_mod1},
    axis=0,
    join="outer",
    label="group",
    fill_value=0,
    index_unique="-"
)
embedder = TruncatedSVD(n_components=par['n_pcs'])
dr_mod1 = embedder.fit_transform(input_train.X)

# split dimred back up
dr_train_mod1 = dr_mod1[input_train.obs['group'] == 'train']
dr_test_mod1 = dr_mod1[input_train.obs['group'] == 'test']


logging.info('train regression model to predict the train_mod2 data from the dr_mod1 values')
reg = KNeighborsRegressor(
    n_neighbors=par['n_neighbors'],
    metric=par['distance_method']
)
reg.fit(dr_train_mod1, input_train_mod2.X.toarray())


logging.info('predict test_mod2 matrix from model and test_mod1')
pred_test_mod2 = reg.predict(dr_test_mod1)


logging.info('calculate k nearest neighbors between test_mod2 and predicted test_mod2')
# To get the matching matrix, for each point in mod1_test, we take the 100 nearest neighbors of that
# point in the transformed mod2_test dataset
n_neighbors = min(100, input_test_mod1.n_obs)
nn = NearestNeighbors(n_neighbors=n_neighbors).fit(pred_test_mod2)
distances, indices = nn.kneighbors(X=input_test_mod2.X)

logging.info('transform k nearest neighbors into a pairing matrix')
# Translate the neighborhood assignments to a pairing matrix that is (n_obs, n_obs)
# NOTE: `pairing_matrix` must have NO MORE than 100*n_obs non-zero entries
ind_i = np.tile(np.arange(input_test_mod1.n_obs), (n_neighbors, 1)).T.flatten()
ind_j = indices.flatten()
ind_dist = distances.flatten()
ind_x = 2 * max(ind_dist) - ind_dist
pairing_matrix = scipy.sparse.csr_matrix(
    (ind_x, (ind_i, ind_j)),
    dims=(input_test_mod1.n_obs, input_test_mod2.n_obs)
)

logging.info('write prediction output')
out = ad.AnnData(
    X=pairing_matrix,
    uns={
        "dataset_id": input_train_mod1.uns["dataset_id"],
        "method_id": method_id
    }
)
out.write_h5ad(par['output'], compression = "gzip")
