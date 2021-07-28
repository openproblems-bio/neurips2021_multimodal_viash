# adapted from https://github.com/YosefLab/totalVI_reproducibility/blob/master/harmonization/sln_harmo_d1_111_d2_206.ipynb

import csv
import gzip
import os
import sys
import scipy.io
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import torch

import scvi
from scvi.dataset import GeneExpressionDataset, CellMeasurement, AnnDatasetFromAnnData
from scvi.models import VAE, TOTALVI
from scvi.inference import TotalPosterior, TotalTrainer, Posterior, UnsupervisedTrainer

import anndata
import scanpy as sc
import umap

from sklearn.preprocessing import normalize

# Control UMAP numba warnings
import warnings; warnings.simplefilter('ignore')

import hotspot
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import normalize

sys.path.append("../utils/")
from utils import *

# for measurment mixing metric
from scipy.stats import mannwhitneyu

# VIASH START
par = {
    "input_mod1": "/resources_test/joint_embedding/test_resource.solution.h5ad",
    "input_mod2": "/resources_test/joint_embedding/test_resource.prediction.h5ad",
    "output": "/resources_test/joint_embedding/test_resource.scores_totalvi.h5ad",
}
# VIASH END

save_path = "/mnt/ibm_lg/achen/neurips/totalvi/"

## TO BE FILLED
#dataset1 = 
#dataset2 = 

#anndata.uns["metric_id"] = 
#anndata.uns["metric_values"]
#%load_ext autoreload
#%autoreload 2
#%matplotlib inline

overwrite=False

set_seed(123)

N_PCS_SEURAT = 30
N_PCS_SCAN = 100

#os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID" 
#os.environ['CUDA_VISIBLE_DEVICES']='1'

colors = ["#3B7EA1", "#FDB515", "#D9661F", "#859438", "#EE1F60", "#00A598"]
# sc.set_figure_params(figsize=(4, 4))
sns.set(context="notebook", font_scale=1.3, style="ticks")
# sns.set(style="ticks")
sns.set_palette(sns.color_palette(colors))
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['savefig.transparent'] = True
plt.rcParams['figure.figsize'] = (4, 4)

sc.settings._vector_friendly = True
DPI = 300
W_SPACE = 0.18

def load_datasets(dataset1, dataset2):
    """
    dataset1 and dataset2 are paths to h5 objects to be integrated. Returns dataset (list) and norm_full_data.
    """
    anndataset_111 = anndata.read(dataset1)
    anndataset_111 = anndataset_111[anndataset_111.obs["batch_indices"] == 0]
    anndataset_206 = anndata.read(dataset2)
    anndataset_206 = anndataset_206[anndataset_206.obs["batch_indices"] == 1]
    anndataset_206.obs["batch_indices"] -= 1

    keep_pro_111 = [not p.startswith("HTO") for p in anndataset_111.uns["protein_names"]]
    keep_pro_206 = [not (p.startswith("HTO") or p.startswith("ADT_Isotype")) for p in anndataset_206.uns["protein_names"]]

    anndataset_111.obsm["protein_expression"] = anndataset_111.obsm["protein_expression"][:, keep_pro_111]
    anndataset_111.uns["protein_names"] = anndataset_111.uns["protein_names"][keep_pro_111]
    anndataset_206.obsm["protein_expression"] = anndataset_206.obsm["protein_expression"][:, keep_pro_206]
    anndataset_206.uns["protein_names"] = anndataset_206.uns["protein_names"][keep_pro_206]

    hvg_111 = anndataset_111.var["hvg_encode"]
    hvg_206 = anndataset_206.var["hvg_encode"]

    assert (hvg_111 == hvg_206).all()

    dataset_111 = AnnDatasetFromAnnData(ad=anndataset_111[:, hvg_111])
    protein_data_111 = CellMeasurement(
        name="protein_expression",
        data=anndataset_111.obsm["protein_expression"].astype(np.float32),
        columns_attr_name="protein_names",
        columns=anndataset_111.uns["protein_names"],
    )
    dataset_111.initialize_cell_measurement(protein_data_111)
    dataset_111.gene_names = anndataset_111[:, hvg_111].var_names.values

    dataset_206 = AnnDatasetFromAnnData(ad=anndataset_206[:, hvg_111])
    protein_data_206 = CellMeasurement(
        name="protein_expression",
        data=anndataset_206.obsm["protein_expression"].astype(np.float32),
        columns_attr_name="protein_names",
        columns=anndataset_206.uns["protein_names"],
    )
    dataset_206.initialize_cell_measurement(protein_data_206)
    dataset_206.gene_names = anndataset_206[:, hvg_206].var_names.values

    dataset = GeneExpressionDataset()
    dataset.populate_from_datasets([dataset_111, dataset_206])

    # standard log library size normalize RNA and concatenate log transformed protein
    full_data = dataset.X
    norm_full_data = anndata.AnnData(X=full_data)
    sc.pp.normalize_per_cell(norm_full_data, counts_per_cell_after=1e4)
    sc.pp.log1p(norm_full_data)
    norm_full_data = anndata.AnnData(
        X=np.concatenate(
            [
                norm_full_data.X,
                np.log1p(dataset.protein_expression),
            ],
            axis=1,
        )
    )
    col_names = np.concatenate([dataset.gene_names, dataset.protein_names])
    norm_full_data.var_names = col_names
#   return dataset_111, dataset_206, dataset
    datasets = [dataset_111, dataset_206, dataset]
    return datasets, norm_full_data

def run_models(datasets):
    """
    Train models on datasets. 
    """
    models = []
    trainers = []
    names = ["SLN_D1_111", "SLN_D2_206", "SLN_harmo_intersect"]

    for d in datasets:
        if d.n_batches > 1:
            m = TOTALVI(d.nb_genes, d.protein_expression.shape[1], n_latent=20, n_batch=d.n_batches, encoder_batch=True)
        else:
            m = TOTALVI(d.nb_genes, d.protein_expression.shape[1], n_latent=20)
        models.append(m)
        use_cuda = True
        lr = 4e-3
        early_stopping_kwargs = {
            "early_stopping_metric": "elbo",
            "save_best_state_metric": "elbo",
            "patience": 45,
            "threshold": 0,
            "reduce_lr_on_plateau": True,
            "lr_patience": 30,
            "lr_factor": 0.6,
            "posterior_class": TotalPosterior,
        }
        
        trainer = TotalTrainer(
            m,
            d,
            train_size=0.90,
            test_size=0.10,
            use_cuda=use_cuda,
            frequency=1,
            data_loader_kwargs={"batch_size":256, "pin_memory":False},
            early_stopping_kwargs=early_stopping_kwargs,
        )
        trainers.append(trainer)
    totalvae_111, totalvae_206, totalvae = models[0], models[1], models[2]

    if overwrite is True:
        for t, n in zip(trainers, names):
            t.train(lr=lr, n_epochs=500)
            torch.save(t.model.state_dict(), "saved_models/" + n + ".pt")
    else:
        for m, n, t in zip(models, names, trainers):
            try:
                m.load_state_dict(torch.load("saved_models/" + n + ".pt"))
                m.eval()
            except FileNotFoundError:
                t.train(lr=lr, n_epochs=500)
                torch.save(t.model.state_dict(), "saved_models/" + n + ".pt")
    # create posterior on full data
    posteriors = []
    latents = []
    for t, m, d in zip(trainers, models, datasets):
        full_posterior = t.create_posterior(
            m, d, indices=np.arange(len(d)), type_class=TotalPosterior
        )
        posteriors.append(full_posterior)
        # extract latent space
        latent_mean, batch_index, label, library_gene = full_posterior.sequential().get_latent()
        latents.append(latent_mean)

    post_adatas = []
    if overwrite:
        for n in names:
            os.remove("saved_post_adata/" + n + ".h5ad")
    for n in names:
        try:
                post_adatas.append(anndata.read_h5ad("saved_post_adata/" + n + ".h5ad"))
        except OSError:
            for d, l, n in zip(datasets, latents, names):
                post_adata = anndata.AnnData(X=d.X)
                post_adata.var.index = d.gene_names
                post_adata.obsm["X_totalVI"] = l
                sc.pp.neighbors(post_adata, use_rep="X_totalVI", n_neighbors=25, metric="correlation")
                sc.tl.umap(post_adata, min_dist=0.2)
                sc.tl.leiden(post_adata, key_added="leiden_totalVI", resolution=0.8)
                post_adata.write("saved_post_adata/" + n + ".h5ad", compression="gzip")
                post_adatas.append(post_adata)
    post_adatas[2].obs["batch_indices"] = [str(b[0]) for b in datasets[2].batch_indices]
    return post_adatas, posteriors

def calculate_metrics(post_adatas, datasets, norm_full_data, posteriors):
    """
    Calculates 4 totalVI metrics on datasets.
    """
    dataset = datasets[2]
    
    ENTROPY_K = 100
    harmo_metrics = pd.DataFrame(
        index=["totalVI-intersect"],
        columns=["Latent mixing metric", "Feature retention metric", "Clustering metric", 
                 'Measurement mixing metric - gene', 'Measurement mixing metric - protein'],
    )
    harmo_metrics["models"] = ["totalVI-intersect"]

    knn = np.ceil(len(dataset.X) * (np.arange(0.1, 1, 0.1)) / 100).astype(int)

    harmo_metrics.loc["totalVI-intersect", "Latent mixing metric"] = entropy_batch_mixing(
        post_adatas[2].obsm["X_totalVI"], datasets[2].batch_indices.ravel(), n_neighbors=ENTROPY_K
    )

    harmo_metrics.loc["totalVI-intersect", "Feature retention metric"], hs_total_1_joint, hs_total_1 = hotspot_score(
        norm_full_data,
        post_adatas[0].obsm["X_totalVI"],
        post_adatas[1].obsm["X_totalVI"],
        post_adatas[2].obsm["X_totalVI"],
        datasets[2].batch_indices.ravel(),
    )

    harmo_metrics.loc["totalVI-intersect", "Clustering metric"] = np.mean(clustering_metric_silhoutte(
        post_adatas[0],
        post_adatas[1],
        post_adatas[2],
        datasets[2].batch_indices.ravel(),
        k=15,
        use_rep="X_totalVI",
        resolution=1.0,
    )
    )

    #measurement mixing metric
    denoised_genes, denoised_proteins = posteriors[-1].sequential().get_normalized_denoised_expression(
        n_samples=50, give_mean=True, transform_batch=[0, 1], 
    )

    denoised_total = np.concatenate([denoised_genes, denoised_proteins], axis=-1)

    totalVI_Us = []
    totalVI_p_vals = []
    mnn_Us = []
    mnn_p_vals = []
    types = []
    models = []
    full_Us = []
    full_ps = []
    mannwhit = pd.DataFrame(columns=["Feature", "Mixing metric", "ps", "Method"])
    batches = datasets[-1].batch_indices.ravel()

    for i in range(denoised_total.shape[1]):
        if i < datasets[2].nb_genes:
            types.append("Gene")
        else:
            types.append("Protein")
        U, p_val = mannwhitneyu(
            denoised_total[batches == 0, i], denoised_total[batches == 1, i]
        )
        totalVI_Us.append(U)
        totalVI_p_vals.append(p_val)
        full_Us.append(U)
        full_ps.append(p_val)
        models.append("totalVI-intersect")
        
    mannwhit["Feature"] = types
    mannwhit["Mixing Metric"] = full_Us
    mannwhit["ps"] = full_ps
    mannwhit["Method"] = models

    #median of measurement mixing metric - gene
    mannwhit_gene_df = mannwhit.loc[mannwhit['Feature'] == 'Gene']
    mixing_metric_gene = np.median(mannwhit_gene_df['Mixing Metric'])

    #median of measurement mixing metric - protein
    mannwhit_protein_df = mannwhit.loc[mannwhit['Feature'] == 'Protein']
    mixing_metric_protein = np.median(mannwhit_protein_df['Mixing Metric'])

    harmo_metrics['Measurement mixing metric - gene'] = mixing_metric_gene
    harmo_metrics['Measurement mixing metric - protein'] = mixing_metric_protein

    return harmo_metrics

def run_totalvi_compute_metrics(dataset1, dataset2):
    """
    Combines functions (combine datasets, run model, and calcuate metrics). Returns totalVI metrics.
    """
    datasets, norm_full_data = load_datasets(dataset1, dataset2)
    post_adatas, posteriors = run_models(datasets)
    harmo_metrics = calculate_metrics(post_adatas, datasets, norm_full_data, posteriors)
    return harmo_metrics

#from totalvi github respository 
#dataset1 = "totalVI_reproducibility/data/spleen_lymph_111.h5ad"
#dataset2 = "totalVI_reproducibility/data/spleen_lymph_206.h5ad"

harmo_metrics = run_totalvi_compute_metrics(dataset1, dataset2)

output = anndata.AnnData(
    X = None, 
    uns={
        "totalvi_metrics": harmo_metrics,
    }
)

output.write_h5ad(par["output"])