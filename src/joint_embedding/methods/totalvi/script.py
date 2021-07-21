print("Load dependencies")

import pandas as pd
import numpy as np
import anndata

import scvi
import scanpy as sc

## VIASH START
par = {
    "input_mod1": "task2/spleen_lymph_111.mod1.h5ad",
    "input_mod2": "task2/spleen_lymph_111.mod2.h5ad",
    "output_prediction": "task2/spleen_lymph_111.mod1.h5ad",
    "hvg_number": 4000,
    "max_epochs": 400,
    
}
## VIASH END

###############################################################################
###                     Load and prepare data                               ###
###############################################################################

print("Load and prepare data")

adata_mod1 = sc.read_h5ad(par['input_mod1'])
adata_mod2 = sc.read_h5ad(par['input_mod2'])


adata_mod1.obsm['protein_expression'] = adata_mod2.X


print('Select highly variable genes')

sc.pp.highly_variable_genes(
    adata_mod1,
    n_top_genes=par['hvg_number'],
    flavor="seurat_v3",
    batch_key="batch",
    subset=True
)


###############################################################################
###                      Set up and train model                             ###
###############################################################################

print("Set up model")

scvi.data.setup_anndata(
    adata_mod1,
    batch_key="batch",
    protein_expression_obsm_key="protein_expression"
)


vae = scvi.model.TOTALVI(adata_mod1, latent_distribution="normal")


print('train totalVI')

vae.train(max_epochs = par['max_epochs'])



###############################################################################
###                             SAVE OUTPUT                                 ###
###############################################################################

print("Postprocessing and saving output")


adata_out = anndata.AnnData(X=vae.get_latent_representation())
uns = { "dataset_id" : adata_mod1.uns["dataset_id"],
        "method_id" : "totalvi"}

adata_out.uns = uns

adata_out.obs = adata_mod1.obs[['batch']]

adata.write_h5ad(par['output_prediction'], compression = "gzip")

