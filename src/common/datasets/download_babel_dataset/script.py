import tarfile
import urllib.request
import tempfile
import anndata as ad
import pandas as pd
import scanpy as sc

## VIASH START

par = {
    "id": "babel_GM12878",
    "input": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE166797&format=file",
    "organism": "human",
    "output_rna": "output_rna.h5ad",
    "output_mod2": "output_mod2.h5ad",
    "rep_n": 1
}

## VIASH END

tempdir = tempfile.mkdtemp()
# with tempfile.TemporaryDirectory() as tempdir:

print("Downloading file from", par['input'])
temptar = f"{tempdir}/GSE166797_RAW.tar"
urllib.request.urlretrieve(par['input'], temptar)

print('Extracting tar file into', temptar)
tf = tarfile.open(temptar)
tf.extractall(tempdir)
tf.close()

print("Reading sample h5")
ad_rep1 = sc.read_10x_h5(tempdir + '/GSM5085810_GM12878_rep1_filtered_feature_bc_matrix.h5', gex_only = False)
ad_rep2 = sc.read_10x_h5(tempdir + '/GSM5085812_GM12878_rep2_filtered_feature_bc_matrix.h5', gex_only = False)

ad_rep1.var_names_make_unique()
ad_rep2.var_names_make_unique()

print("Merging into one AnnData object")
comb = ad.concat(
    {"rep1": ad_rep1, "rep2": ad_rep2},
    axis=0,
    join="outer",
    label="batch",
    fill_value=0,
    index_unique="-"
)
catvar = pd.concat([ad_rep1.var, ad_rep2.var]).drop_duplicates()
comb.var = catvar.loc[comb.var_names]

print("Splitting up into GEX and ATAC")
adata_RNA = comb[:, comb.var.feature_types == "Gene Expression"].copy()
adata_ATAC = comb[:, comb.var.feature_types == "Peaks"].copy()

adata_RNA.var.feature_types = 'GEX'
adata_ATAC.var.feature_types = 'ATAC'

uns = { "dataset_id" : par["id"], "organism" : par["organism"] }
adata_RNA.uns = uns
adata_ATAC.uns = uns

print("Writing output")
adata_RNA.write_h5ad(par['output_rna'], compression = "gzip")
adata_ATAC.write_h5ad(par['output_mod2'], compression = "gzip")
