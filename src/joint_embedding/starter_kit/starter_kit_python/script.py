# Dependencies:
# pip: scikit-learn, anndata, scanpy, umap-learn
#
# Python starter kit for the NeurIPS 2021 Single-Cell Competition. Parts
# with `TODO` are supposed to be changed by you.
#
# More documentation:
#
# https://viash.io/docs/creating_components/python/

import logging
import anndata as ad
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
    'n_dim': 100,
}

## VIASH END

# TODO: change this to the name of your method
method_id = "python_starter_kit"

logging.info('Reading `h5ad` files...')

ad_mod1 = ad.read_h5ad(par['input_mod1'])
ad_mod2 = ad.read_h5ad(par['input_mod2'])

# TODO: implement own method

logging.info('Concatenating modality 1 and modality 2...')

# Transpose the objects (concatenate only works along axis 0)
ad_mod1.transpose()
ad_mod2.transpose()

# Concatenate the objects
ad_combined = ad_mod1.concatenate(ad_mod2)

# Transpose back
ad_combined.transpose()

logging.info('Performing dimensionality reduction on concatenated datasets...')

embedder = umap.UMAP(
    n_components=par['n_dim'],
)

X_umap = embedder.fit_transform(ad_combined.X)

# Replace ad.X with the dimensionaltiy reduced embedding
adata = ad.AnnData(
    X=X_umap,
    uns={
        'dataset_id': ad_mod1.uns['dataset_id'],
        'method_id': method_id,
    },
)

logging.info('Storing annotated data...')

adata.write(par['output'])
