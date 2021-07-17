print("Load dependencies")
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
    "id": "azimuth_ref",
    "input_count": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file",
    "input_meta": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file&file=GSE164378%5Fsc%2Emeta%2Edata%5F3P%2Ecsv%2Egz",
    "output_rna": "output_rna.h5ad",
    "output_mod2": "output_mod2.h5ad"
}
## VIASH END

###############################################################################
###                     DOWNLOAD AND READ DATA.                             ###
###############################################################################
print("Downloading file from", par['input_count'])
tar_temp = tempfile.NamedTemporaryFile()
url = par['input_count']
urllib.request.urlretrieve(url, tar_temp.name)

print("Downloading meta data from", par['input_meta'])
meta_temp = tempfile.NamedTemporaryFile()
url = par['input_meta']
urllib.request.urlretrieve(url, meta_temp.name)

###############################################################################
###                      EXTRACT AND CREATE H5ADs                           ###
###############################################################################

print("Extracting and create h5ads")
samples = ['GSM5008737_RNA_3P', 'GSM5008738_ADT_3P'] # first sample is rna, second is protein data

adatas = []        
with tarfile.open(tar_temp.name) as tar:

    for sample in samples:
        print("Processing sample " + sample)
        with gzip.open(tar.extractfile(sample + '-matrix.mtx.gz'), 'rb') as mm:
            print('Loading matrix')
            X = scipy.io.mmread(mm).T.tocsr()
        obs = pd.read_csv(
            tar.extractfile(sample + '-barcodes.tsv.gz'), 
            compression='gzip',
            header=None, 
            sep='\t',
            index_col=0
        )
        obs.index.name = None
        var = pd.read_csv(
            tar.extractfile(sample + '-features.tsv.gz'), 
            compression='gzip',
            header=None, 
            sep='\t'
        ).iloc[:, :1]
        var.columns = ['names']
        var.index = var['names'].values
        adata = anndata.AnnData(X=X, obs=obs, var=var)

        adata.var_names_make_unique()
        adatas.append(adata)

    tar.close()

adata = adatas[0]
protein = adatas[1]

###############################################################################
###                            POST PROCESS                                 ###
###############################################################################
print("Reading metadata")
meta = pd.read_csv(meta_temp.name, index_col = 0, compression = "gzip")
meta_adt = meta.loc[:,~meta.columns.str.endswith('RNA')]
meta_rna = meta.loc[:,~meta.columns.str.endswith('ADT')]

print("Setting additional output fields")
# set obs
adata.obs = adata.obs.join(meta_rna).rename(columns = {'Batch': 'batch'}, inplace = True)
adata.obs['cell_type'] = adata.obs['celltype.l2']

protein.obs = protein.obs.join(meta_adt).rename(columns = {'Batch': 'batch'}, inplace = True)
protein.obs['cell_type'] = protein.obs['celltype.l2']

#  set var
adata.var['feature_types'] = "GEX"
protein.var['feature_types'] = "ADT"

# set uns 
uns = { "dataset_id" : par["id"] }
adata.uns = uns
protein.uns = uns


###############################################################################
###                             SAVE OUTPUT                                 ###
###############################################################################
print("Saving output")
adata.write_h5ad(par['output_rna'], compression = "gzip")
protein.write_h5ad(par['output_mod2'], compression = "gzip")

