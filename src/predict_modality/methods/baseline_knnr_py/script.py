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

from scipy.sparse import csc_matrix

from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import KNeighborsRegressor

logging.basicConfig(level=logging.INFO)

## VIASH START
# Anything within this block will be removed by `viash` and will be
# replaced with the parameters as specified in your config.vsh.yaml.
par = {
    'input_train_mod1': 'sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod1.h5ad',
    'input_train_mod2': 'sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod2.h5ad',
    'input_test_mod1': 'sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod1.h5ad',
    'distance_method': 'minkowski',
    'output': 'output.h5ad',
    'n_pcs': 4,
    'n_neighbors': 5,
}
meta = { 'functionality_name': 'foo' }
## VIASH END

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

logging.info('Performing dimensionality reduction on modality 1 values...')
embedder = TruncatedSVD(n_components=par['n_pcs'])
X = embedder.fit_transform(input_train.X)

# split dimred back up
X_train = X[input_train.obs['group'] == 'train']
X_test = X[input_train.obs['group'] == 'test']
y_train = input_train_mod2.X.toarray()

assert len(X_train) + len(X_test) == len(X)

logging.info('Running KNN regression...')

reg = KNeighborsRegressor(
    n_neighbors=par['n_neighbors'],
    metric=par['distance_method']
)

reg.fit(X_train, y_train)
y_pred = reg.predict(X_test)

y_pred = csc_matrix(y_pred)

adata = ad.AnnData(
    X=y_pred,
    obs=input_test_mod1.obs,
    var=input_train_mod2.var,
    uns={
        'dataset_id': input_train_mod1.uns['dataset_id'],
        'method_id': meta["functionality_name"],
    },
)

logging.info('Storing annotated data...')
adata.write_h5ad(par['output'], compression = "gzip")
