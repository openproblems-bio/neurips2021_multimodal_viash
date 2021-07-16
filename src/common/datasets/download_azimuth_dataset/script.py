

## VIASH START
par = {
    "id": "azimuth_ref",
    "input_count": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file",
    "input_meta": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file&file=GSE164378%5Fsc%2Emeta%2Edata%5F3P%2Ecsv%2Egz",
    "out_raw_count": "resources_test_2/common/azimuth_raw/GSE164378_RAW.tar",
    "out_raw_meta": "resources_test_2/common/azimuth_raw/GSE164378_sc.meta.data_3P.csv.gz",
    "output_rna": "output_rna.h5ad",
    "output_mod2": "output_mod2.h5ad"
}
## VIASH END



###############################################################################
###                            LOAD DEPENDENCIES                            ###
###############################################################################

import urllib.request
import tempfile
import anndata
import scanpy as sc
import pandas as pd

import tarfile
import numpy as np
import gzip
import scipy.io

###############################################################################
###                         DOWNLAOD RAW FILES                              ###
###############################################################################


print("Downloading count matrices from", par['input_count'])

url = par['input_count']
urllib.request.urlretrieve(url, par['out_raw_count'])


print("Downloading meta data from", par['input_meta'])

url = par['input_meta']
urllib.request.urlretrieve(url, par['out_raw_meta'])


###############################################################################
###                      EXTRACT AND CREATE H5ADs                           ###
###############################################################################

fn = par['out_raw_count']

adatas = []        
with tarfile.open(fn) as tar:
    samples = ['GSM5008737_RNA_3P', 'GSM5008738_ADT_3P'] # first sample is rna, second is protein data

    for sample in samples:
        print(sample)
        with gzip.open(tar.extractfile(sample + '-matrix.mtx.gz'), 'rb') as mm:
            print('Loading matrix')
            X = scipy.io.mmread(mm).T.tocsr()
        obs = pd.read_csv(tar.extractfile(sample + '-barcodes.tsv.gz'), compression='gzip',
                        header=None, sep='\t', index_col=0)
        obs.index.name = None
        var = pd.read_csv(tar.extractfile(sample + '-features.tsv.gz'), compression='gzip',
                        header=None, sep='\t').iloc[:, :1]
        var.columns = ['names']
        var.index = var['names'].values
        adata = anndata.AnnData(X=X, obs=obs, var=var)

        adata.var_names_make_unique()
        adatas.append(adata)

    tar.close()

adata = adatas[0]
protein = adatas[1]



meta = pd.read_csv(par['out_raw_meta'], index_col =0)

meta_adt = meta.loc[:,~meta.columns.str.endswith('RNA')]

meta_rna = meta.loc[:,~meta.columns.str.endswith('ADT')]

adata.obs = protein.obs.join(meta_rna)
protein.obs = protein.obs.join(meta_adt)

adata.var['feature_types'] = "GEX"
protein.var['feature_types'] = "ADT"

adata.uns = {'modality':'GEX'}
protein.uns = {'modality':'ADT'}



###############################################################################
###                             SAVE OUTPUT                                 ###
###############################################################################


adata.write_h5ad(par['output_rna'], compression = "gzip")
protein.write_h5ad(par['output_mod2'], compression = "gzip")
