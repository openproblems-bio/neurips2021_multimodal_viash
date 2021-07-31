import scvi
from scvi.dataset import GeneExpressionDataset, CellMeasurement, AnnDatasetFromAnnData
from scvi.models import VAE, TOTALVI
from scvi.inference import TotalPosterior, TotalTrainer, Posterior, UnsupervisedTrainer
import os
import sys

save_path = "/mnt/ibm_lg/achen/neurips/totalvi/"

import anndata
import scanpy as sc
import numpy as np
import pandas as pd

sys.path.append("../utils/")
from utils import *

# VIASH START
par = {
    "input_prediction": "/resources_test/joint_embedding/test_resource.prediction.h5ad",
    "input_solution": "/resources_test/joint_embedding/test_resource.solution.h5ad",
    "output": "/resources_test/joint_embedding/test_resource.scores_totalvi.h5ad",
}
# VIASH END

predict_adata = sc.read_h5ad("input_prediction")
solution_mod2 = sc.read_h5ad("input_solution")

method_name = predict_adata.uns['method_id'] 

metrics_adata = predict_adata
metrics_adata.obsm[method_name] = predict_adata.X

ENTROPY_K = 100
harmo_metrics = pd.DataFrame(
    index=[method_name],
    columns=["Latent mixing metric", "Feature retention metric", "Clustering metric", 
             'Measurement mixing metric - gene', 'Measurement mixing metric - protein'],

knn = np.ceil(metrics_adata.n_obs * (np.arange(0.1, 1, 0.1)) / 100).astype(int)

# calculate latent mixing metric
latent_mixing_metric = entropy_batch_mixing(
    metrics_adata.obsm[method_name], metrics_adata.obs["batch"].values, n_neighbors=ENTROPY_K
)

#save output
adata_out = anndata.AnnData()

uns = {"dataset_id": predict_adata.uns["dataset_id"],
      "method_id" : predict_adata.uns["method_id"],
      "method_ids" : ["latent mixing metric"],
      "metric_values" : [latent_mixing_metric],
      "metric_moreisbetter" : [True]}

adata_out.uns = uns

adata_out.write_h5ad('output_metrics', compression = "gzip")