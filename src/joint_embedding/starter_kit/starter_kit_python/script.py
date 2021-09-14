# Dependencies:
# pip: anndata, umap-learn
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

from sklearn.decomposition import TruncatedSVD

logging.basicConfig(level=logging.INFO)

## VIASH START

# Anything within this block will be removed by `viash` and will be
# replaced with the parameters as specified in your config.vsh.yaml.

dataset_path = 'sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.'
# dataset_path = 'output/datasets/joint_embedding/openproblems_bmmc_multiome_phase1/openproblems_bmmc_multiome_phase1.censor_dataset.output_'

par = {
    'input_mod1': dataset_path + 'mod1.h5ad',
    'input_mod2': dataset_path + 'mod2.h5ad',
    'output': 'output.h5ad',
    'n_dim': 50,
}

## VIASH END

# TODO: change this to the name of your method
method_id = "python_starter_kit"

logging.info('Reading `h5ad` files...')
ad_mod1 = ad.read_h5ad(par['input_mod1'])
ad_mod2 = ad.read_h5ad(par['input_mod2'])

# TODO: implement your own method
logging.info('Performing dimensionality reduction on modality 1 values...')
embedder_mod1 = TruncatedSVD(n_components=int(par["n_dim"]/2))
mod1_pca = embedder_mod1.fit_transform(ad_mod1.X)
mod1_obs = ad_mod1.obs
mod1_uns = ad_mod1.uns
del ad_mod1

logging.info('Performing dimensionality reduction on modality 2 values...')
embedder_mod1 = TruncatedSVD(n_components=int(par["n_dim"]/2))
mod2_pca = embedder_mod1.fit_transform(ad_mod2.X)
del ad_mod2

logging.info('Concatenating datasets')
pca_combined = np.concatenate([mod1_pca, mod2_pca], axis=1)

logging.info('Storing output to file')
adata = ad.AnnData(
    X=pca_combined,
    obs=mod1_obs,
    uns={
        'dataset_id': mod1_uns['dataset_id'],
        'method_id': method_id,
    },
)
adata.write_h5ad(par['output'], compression="gzip")
