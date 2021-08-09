import anndata as ad

# VIASH START
par = {
    "input_prediction": "resources_test/joint_embedding/test_resource.prediction.h5ad",
    "input_solution": "resources_test/joint_embedding/test_resource.solution.h5ad",
    "output": "resources_test/joint_embedding/test_resource.scores_totalvi.h5ad",
    "n_neighbors": 100
}
resources_dir = "src/joint_embedding/metrics/calculate_totalVI_metrics/"
# VIASH END

print("Import utils library")
import sys
sys.path.append(resources_dir)
from utils import entropy_batch_mixing

print("Read input files")
predict_adata = ad.read_h5ad(par["input_prediction"])
solution_adata = ad.read_h5ad(par["input_solution"])

print("Merge prediction with solution")
merged_adata = predict_adata.copy()
merged_adata.obs["batch"] = solution_adata.obs["batch"][merged_adata.obs_names]

print("Calculate latent mixing metric")
latent_mixing_metric = entropy_batch_mixing(
    merged_adata.X,
    merged_adata.obs["batch"].values,
    n_neighbors=par["n_neighbors"]
)

print("Write output")
adata_out = ad.AnnData(
    uns = {
        "dataset_id": predict_adata.uns["dataset_id"],
        "method_id" : predict_adata.uns["method_id"],
        "method_ids" : ["latent_mixing_metric"],
        "metric_values" : [latent_mixing_metric]
    }
)

adata_out.write_h5ad(par['output'], compression = "gzip")