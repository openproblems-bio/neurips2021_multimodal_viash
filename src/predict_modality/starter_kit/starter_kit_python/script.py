# Dependencies:
# pip: scikit-learn, anndata, scanpy
#
# Python starter kit for the NeurIPS 2021 Single-Cell Competition. Parts
# with `TODO` are supposed to be changed by you.
#
# More documentation:
#
# https://viash.io/docs/creating_components/python/

import anndata
import scanpy
import logging

from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
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
    'metric': 'minkowski',
    'output': 'output.h5ad',
    'n_components': 4,
    'n_neighbors': 5,
}

## VIASH END

# TODO: change this to the name of your method
method_id = "python_starter_kit"

logging.info('Reading `h5ad` files...')

data_modality_1 = scanpy.read_h5ad(par['input_mod1'])
data_modality_2 = scanpy.read_h5ad(par['input_mod2'])

logging.info('Performing dimensionality reduction on modality 1 values...')

mds = MDS(
    n_components=par['n_components'],
    dissimilarity='precomputed',
)

D = pairwise_distances(
    data_modality_1.X.toarray(), metric=par['metric']
)

X = mds.fit_transform(D)

X_train = X[data_modality_1.obs['group'] == 'train']
X_test = X[data_modality_1.obs['group'] == 'test']

assert len(X_train) + len(X_test) == len(X)

# Get all responses of the training data set. If you are
# a `scikit-learn` user the nomenclature may be somewhat
# unknown.
#
# Make sure to use `toarray()` because the output might
# be sparse and `KNeighborsRegressor` cannot handle it.
y_train = data_modality_2.X.toarray()

logging.info('Running KNN regression...')

reg = KNeighborsRegressor(
    n_neighbors=par['n_neighbors'],
    metric=par['metric']
)

reg.fit(X_train, y_train)
y_pred = reg.predict(X_test)

adata = anndata.AnnData(
    X=y_pred,
    uns={
        'dataset_id': data_modality_1.uns['dataset_id'],
        'method_id': method_id,
    },
    # FIXME: store `obs_names` and `var_names` also?
)

logging.info('Storing annotated data...')

adata.write(par['output'])
