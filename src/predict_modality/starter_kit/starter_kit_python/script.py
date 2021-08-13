# Dependencies:
# pip: scikit-learn, anndata, scanpy
#
# Python starter kit for the NeurIPS 2021 Single-Cell Competition. Parts
# with `TODO` are supposed to be changed by you.
#
# More documentation:
#
# https://viash.io/docs/creating_components/python/

import logging
import anndata as ad

from scipy.sparse import csc_matrix

from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import KNeighborsRegressor

logging.basicConfig(level=logging.INFO)

## VIASH START
# Anything within this block will be removed by `viash` and will be
# replaced with the parameters as specified in your config.vsh.yaml.
par = {
    'input_train_mod1': 'sample_data/test_resource.train_mod1.h5ad',
    'input_train_mod2': 'sample_data/test_resource.train_mod2.h5ad',
    'input_test_mod1': 'sample_data/test_resource.test_mod1.h5ad',
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
input_test_mod1 = ad.read_h5ad(par['input_test_mod1'])

input_train = ad.concat(
    {"train": input_train_mod1, "test": input_test_mod1},
    axis=0,
    join="outer",
    label="group",
    fill_value=0,
    index_unique="-"
)

# TODO: implement own method

logging.info('Performing dimensionality reduction on modality 1 values...')
embedder = TruncatedSVD(n_components=par['n_pcs'])
X = embedder.fit_transform(input_train.X)

# split dimred back up
X_train = X[input_train.obs['group'] == 'train']
X_test = X[input_train.obs['group'] == 'test']
y_train = input_train_mod2.X.toarray()

assert len(X_train) + len(X_test) == len(X)

# Get all responses of the training data set to fit the
# KNN regressor later on.
#
# Make sure to use `toarray()` because the output might
# be sparse and `KNeighborsRegressor` cannot handle it.

logging.info('Running KNN regression...')

reg = KNeighborsRegressor(
    n_neighbors=par['n_neighbors'],
    metric=par['distance_method']
)

reg.fit(X_train, y_train)
y_pred = reg.predict(X_test)

# Store as sparse matrix to be efficient. Note that this might require
# different classifiers/embedders before-hand. Not every class is able
# to support such data structures.
y_pred = csc_matrix(y_pred)

adata = ad.AnnData(
    X=y_pred,
    uns={
        'dataset_id': input_train_mod1.uns['dataset_id'],
        'method_id': method_id,
    },
)

logging.info('Storing annotated data...')
adata.write_h5ad(par['output'], compression = "gzip")
