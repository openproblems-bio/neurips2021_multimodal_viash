import scanpy as sc
import anndata
from utils import entropy_batch_mixing

# VIASH START
par = {
    "input_prediction": "/resources_test/joint_embedding/test_resource.prediction.h5ad",
    "input_solution": "/resources_test/joint_embedding/test_resource.solution.h5ad",
    "output": "/resources_test/joint_embedding/test_resource.scores_totalvi.h5ad",
}
# VIASH END

predict_adata = sc.read_h5ad(par["input_prediction"])
solution_adata = sc.read_h5ad(par["input_solution"])

method_name = predict_adata.uns["method_id"] 

metrics_adata = predict_adata
metrics_adata.obsm[method_name] = predict_adata.X

# calculate mixing statistics
ENTROPY_K = 100

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

adata_out.write_h5ad(par['output'], compression = "gzip")