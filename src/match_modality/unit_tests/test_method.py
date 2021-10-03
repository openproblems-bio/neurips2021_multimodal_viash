from os import path
import subprocess
import anndata as ad
import numpy as np
from scipy.sparse import issparse

## VIASH START
# This code block will be replaced by viash at runtime.
meta = { 'functionality_name': 'foo' }
## VIASH END

method_id = meta['functionality_name']
command = "./" + method_id

# define some filenames
testpar = {
    "input_train_mod1": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod1.h5ad",
    "input_train_mod2": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod2.h5ad",
    "input_train_sol": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_sol.h5ad",
    "input_test_mod1": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod1.h5ad",
    "input_test_mod2": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod2.h5ad",
    "input_test_sol": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_sol.h5ad",
    "output": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
}

print("> Running method")
out = subprocess.check_output([
  command,
  "--input_train_mod1", testpar['input_train_mod1'],
  "--input_train_mod2", testpar['input_train_mod2'],
  "--input_train_sol", testpar['input_train_sol'],
  "--input_test_mod1", testpar['input_test_mod1'],
  "--input_test_mod2", testpar['input_test_mod2'],
  "--output", testpar['output']
]).decode("utf-8")

print("> Checking whether output files were created")
assert path.exists(testpar['output'])

print("> Reading h5ad files")
ad_sol = ad.read_h5ad(testpar['input_test_sol'])
ad_pred = ad.read_h5ad(testpar['output'])

print("> Checking dataset id")
assert ad_pred.uns['dataset_id'] == ad_sol.uns['dataset_id']

print("> Checking method id")
assert ad_pred.uns['method_id'] == method_id

print("> Checking X")
assert issparse(ad_pred.X)
assert np.all([x >= 0 for x in ad_pred.X.nonzero()]), "Values must be strictly non-negative."
assert ad_pred.X.nonzero()[0].size <= 1000 * ad_sol.n_obs
assert ad_pred.n_obs == ad_sol.n_obs
assert ad_pred.n_vars == ad_sol.n_vars
assert np.isclose(ad_pred.X.sum(axis=1), 1, atol=1e-10).all(), "All rows should sum to 1."

print("> Test succeeded!")