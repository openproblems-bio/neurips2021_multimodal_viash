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
import umap.umap_ as umap

logging.basicConfig(level=logging.INFO)

## VIASH START

# Anything within this block will be removed by `viash` and will be
# replaced with the parameters as specified in your config.vsh.yaml.

par = {
    'input_mod1': 'sample_data/test_resource.mod1.h5ad',
    'input_mod2': 'sample_data/test_resource.mod2.h5ad',
    'output': 'output.h5ad',
    'n_dim': 100,
}

## VIASH END

# TODO: change this to the name of your method
method_id = "python_starter_kit"

logging.info('Reading `h5ad` files...')
ad_mod1 = ad.read_h5ad(par['input_mod1'])
ad_mod2 = ad.read_h5ad(par['input_mod2'])

# TODO: implement your own method
logging.info('Concatenating modality 1 and modality 2')
ad_combined = ad.concat([ad_mod1, ad_mod2], axis=1)

logging.info('Performing dimensionality reduction on concatenated datasets')
embedder = umap.UMAP(n_components=par['n_dim'])
X_umap = embedder.fit_transform(ad_combined.X)

logging.info('Storing output to file')
adata = ad.AnnData(
    X=X_umap,
    obs=ad_mod1.obs,
    uns={
        'dataset_id': ad_mod1.uns['dataset_id'],
        'method_id': method_id,
    },
)
adata.write_h5ad(par['output'], compression="gzip")
