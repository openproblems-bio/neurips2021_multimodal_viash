# Dependencies:
# pip: scikit-learn, anndata, scanpy, umap-learn
#
# Python starter kit for the NeurIPS 2021 Single-Cell Competition. Parts
# with `TODO` are supposed to be changed by you.
#
# More documentation:
#
# https://viash.io/docs/creating_components/python/

import anndata
import logging
import scanpy
import umap

from scipy.sparse import csc_matrix

from sklearn.neighbors import KNeighborsRegressor


logging.basicConfig(
    level=logging.INFO
)


## VIASH START

# Anything within this block will be removed by `viash` and will be
# replaced with the parameters as specified in your config.vsh.yaml.

par = {
    'input_mod1': 'sample_data/test_resource.mod1.h5ad',
    'input_mod2': 'sample_data/test_resource.mod2.h5ad',
    'distance_method': 'minkowski',
    'output': 'output.h5ad',
    'n_pcs': 4,
    'n_neighbors': 5,
}

## VIASH END

# TODO: change this to the name of your method
method_id = "python_starter_kit"

logging.info('Reading `h5ad` files...')

data_modality_1 = scanpy.read_h5ad(par['input_mod1'])
data_modality_2 = scanpy.read_h5ad(par['input_mod2'])

logging.info('Performing dimensionality reduction on modality 1 values...')

# Notice how this instantiation also uses the pre-defined parameter for
# the distance method to be used here.
embedder = umap.UMAP(
    n_components=par['n_pcs'],
    metric=par['distance_method'],
)

X = embedder.fit_transform(data_modality_1.X)

X_train = X[data_modality_1.obs['group'] == 'train']
X_test = X[data_modality_1.obs['group'] == 'test']

assert len(X_train) + len(X_test) == len(X)

# Get all responses of the training data set to fit the
# KNN regressor later on.
#
# Make sure to use `toarray()` because the output might
# be sparse and `KNeighborsRegressor` cannot handle it.
y_train = data_modality_2.X.toarray()

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

adata = anndata.AnnData(
    X=y_pred,
    uns={
        'dataset_id': data_modality_1.uns['dataset_id'],
        'method_id': method_id,
    },
)

logging.info('Storing annotated data...')

adata.write(par['output'])
