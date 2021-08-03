import urllib.request
import tempfile
import anndata
import scanpy as sc
import pandas as pd

import tarfile
import numpy as np
import gzip
import scipy.io


## VIASH START

par = {
    "id": "babel_GM12878",
    "input_count": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE166797&format=file",
    "output_rna": "output_rna.h5ad",
    "output_mod2": "output_mod2.h5ad"
    "rep_n": 1
}

## VIASH END


print("Downloading file from", par['input_count'])
tar_temp = tempfile.NamedTemporaryFile()
url = par['input_count']

urllib.request.urlretrieve(url, tar_temp.name)

# two replicates available as independent datasets
samples = ['GSM5085810_GM12878_rep1',
          'GSM5085812_GM12878_rep2']

sample = samples[par['rep_n']]

with tarfile.open(tar_temp.name) as tar:

    with tempfile.TemporaryDirectory() as tmpdirname:
        print('Extracting tar file into:', tmpdirname)
        print("Processing sample " + sample)
        # put tmp directory here instead:
        tar.extractall(path=tmpdirname+'/')

        #adata_file = tar.extractfile('tmp/',sample + '_filtered_feature_bc_matrix.h5')
        adata = sc.read_10x_h5(tmpdirname+'/'+sample + '_filtered_feature_bc_matrix.h5', gex_only = False)


    tar.close()


# Set up var data
adata.var_names_make_unique()

adata.var.feature_types = adata.var.feature_types.astype(str)

adata_RNA = adata[:, adata.var.feature_types.str.startswith('Gene')].copy()

adata_ATAC = adata[:, adata.var.feature_types.str.startswith('Peaks')].copy()

adata_RNA.var.feature_types = 'GEX'
adata_ATAC.var.feature_types = 'ATAC'

uns = { "dataset_id" : par["id"] }
adata_RNA.uns = uns
adata_ATAC.uns =uns

print("Saving output")
adata_RNA.write_h5ad(par['output_rna'], compression = "gzip")
adata_RNA.write_h5ad(par['output_mod2'], compression = "gzip")
