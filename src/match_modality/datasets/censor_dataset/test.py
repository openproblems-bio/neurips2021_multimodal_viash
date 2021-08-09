import subprocess
from os import path

import anndata as ad

print("> Running censor component")
cmd_pars = [
  "./censor_dataset",
  "--input", "dataset.h5ad", 
  "--output_censored", "output_censored.h5ad",
  "--output_solution", "output_solution.h5ad"
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

print("> Checking whether output files were created")
assert path.exists("output_censored.h5ad"), "No censored output was created."
assert path.exists("output_solution.h5ad"), "No solution output was created."

print("> Checking contents of output_censored.h5ad")
adata_orig = ad.read_h5ad("dataset.h5ad")
adata_censor = ad.read_h5ad("output_censored.h5ad")

assert adata_censor.uns["dataset_id"] == (adata_orig.uns["dataset_id"] + "_task2")
assert adata_censor.n_obs == adata_orig.n_obs
assert adata_censor.n_vars == adata_orig.n_vars
assert "modality1" in adata_censor.layers.keys()
assert "modality2" in adata_censor.layers.keys()

print("> Checking contents of output_solution.h5ad")
adata_sol = ad.read_h5ad("output_solution.h5ad")
assert adata_sol.uns["dataset_id"] == (adata_orig.uns["dataset_id"] + "_task2")
assert adata_sol.n_obs == adata_orig.n_obs
assert adata_sol.n_vars == adata_orig.n_vars
assert "shuffle_index" in adata_sol.obs_keys()

# TODO: add more checks (check whether shuffle indices are correct)

print("> Test succeeded!")
