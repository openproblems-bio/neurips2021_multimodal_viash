from os import path
import subprocess
import anndata as ad
from scipy.sparse import issparse

par = {
    "input_train_mod1": "resources_test/match_modality/test_resource.train_mod1.h5ad",
    "input_train_mod2": "resources_test/match_modality/test_resource.train_mod2.h5ad",
    "input_train_sol": "resources_test/match_modality/test_resource.train_sol.h5ad",
    "input_test_mod1": "resources_test/match_modality/test_resource.test_mod1.h5ad",
    "input_test_mod2": "resources_test/match_modality/test_resource.test_mod2.h5ad",
    "input_test_sol": "resources_test/match_modality/test_resource.test_sol.h5ad",
    "output": "output.h5ad",
}

print("> Running method")
out = subprocess.check_output([
  "./" + meta['functionality_name'],
  "--input_train_mod1", par['input_train_mod1'],
  "--input_train_mod2", par['input_train_mod2'],
  "--input_train_sol", par['input_train_sol'],
  "--input_test_mod1", par['input_test_mod1'],
  "--input_test_mod2", par['input_test_mod2'],
  "--output", par['output']
]).decode("utf-8")

print("> Checking whether output files were created")
assert path.exists(par['output'])

print("> Reading h5ad files")
ad_sol = ad.read_h5ad(par['input_test_sol'])
ad_pred = ad.read_h5ad(par['output'])

print("> Checking dataset id")
assert ad_pred.uns['dataset_id'] == ad_sol.uns['dataset_id']

print("> Checking method id")
assert ad_pred.uns['method_id'] == meta['functionality_name']

print("> Checking X")
assert issparse(ad_pred.X)
assert ad_pred.n_obs == ad_sol.n_obs
assert ad_pred.n_vars == ad_sol.n_vars
assert all(ad_pred.obs_names == ad_sol.obs_names)
assert all(ad_pred.var_names == ad_sol.var_names)

print("> Test succeeded!")