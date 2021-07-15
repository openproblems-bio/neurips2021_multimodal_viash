

## VIASH START
par = {
    "id": "totalVI_10x_malt_10k",
    "input": "https://github.com/YosefLab/totalVI_reproducibility/raw/master/data/malt_10k_protein_v3.h5ad",
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

###############################################################################
###                     DOWNLAOD AND READ DATA.                             ###
###############################################################################


print("Downloading file from", par['input'])

h5ad_temp = tempfile.NamedTemporaryFile()
url = par['input']
urllib.request.urlretrieve(url, h5ad_temp.name)


print("Reading h5ad file")

adata = sc.read_h5ad(h5ad_temp.name)

h5ad_temp.close()


###############################################################################
###                     CREATE H5AD FOR BOTH MODALITIES                     ###
###############################################################################


print("Creating h5ad file per modality")


# Extract RNA counts

adata_rna = anndata.AnnData(X = adata.X,
                            obs = adata.obs.loc[:,['n_genes', 'percent_mito', 'n_counts']],
                            var = adata.var,
                            uns = {'modality':'GEX'},
                           )

adata_rna.var['feature_types'] = "GEX"

# Extract ADT (antibody derived transcripts) counts

adata_adt = anndata.AnnData(X = adata.obsm['protein_expression'],
                            var = pd.DataFrame(index=list(adata.uns['protein_names'])),
                            uns = {'modality':'ADT'},
                           )

adata_adt.obs.index = adata.obs.index
adata_adt.var['feature_types'] = "ADT"


###############################################################################
###                             SAVE OUTPUT                                 ###
###############################################################################


adata_rna.write_h5ad(par['output_rna'], compression = "gzip")
adata_adt.write_h5ad(par['output_mod2'], compression = "gzip")


