import anndata
import scvi
import scanpy as sc

## VIASH START
par = {
    "input_mod1": "output/joint_embedding/totalvi_spleen_lymph_111/totalvi_spleen_lymph_111.censor_dataset.output_mod1.h5ad",
    "input_mod2": "output/joint_embedding/totalvi_spleen_lymph_111/totalvi_spleen_lymph_111.censor_dataset.output_mod2.h5ad",
    "output": "tmp/output_prediction.h5ad",
    "hvg_number": 4000,
    "max_epochs": 20
}
## VIASH END

print("Load and prepare data")
adata_mod1 = anndata.read_h5ad(par['input_mod1'])
adata_mod2 = anndata.read_h5ad(par['input_mod2'])
adata_mod1.obsm['protein_expression'] = adata_mod2.X.toarray()

print('Select highly variable genes')
sc.pp.highly_variable_genes(
    adata_mod1,
    n_top_genes=par['hvg_number'],
    flavor="seurat_v3",
    batch_key="batch",
    subset=True
)

print("Set up model")
scvi.data.setup_anndata(
    adata_mod1,
    batch_key="batch",
    protein_expression_obsm_key="protein_expression"
)

print('Train totalVI with', par['max_epochs'], 'epochs')
vae = scvi.model.TOTALVI(adata_mod1, latent_distribution="normal")
vae.train(max_epochs = par['max_epochs'])

print("Postprocessing and saving output")
adata_out = anndata.AnnData(
    X=vae.get_latent_representation(),
    obs=adata_mod1.obs[['batch']],
    uns={
        "dataset_id": adata_mod1.uns["dataset_id"],
        "method_id": "totalvi"
    }
)
adata_out.write_h5ad(par['output'], compression = "gzip")
