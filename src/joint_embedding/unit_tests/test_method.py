from os import path
import subprocess
import anndata as ad

## VIASH START
# This code block will be replaced by viash at runtime.
meta = { 'functionality_name': 'foo' }
## VIASH END

method_id = meta['functionality_name']
command = "./" + method_id

# define some filenames
testpar = {
    "input_mod1": "resources_test/joint_embedding/test_resource.mod1.h5ad",
    "input_mod2": "resources_test/joint_embedding/test_resource.mod2.h5ad",
    "output": "output.h5ad",
}

print("> Running method")
out = subprocess.check_output([
  command,
  "--input_mod1", testpar['input_mod1'],
  "--input_mod2", testpar['input_mod2'],
  "--output", testpar['output']
]).decode("utf-8")

print("> Checking whether output files were created")
assert path.exists(testpar['output'])

print("> Reading h5ad files")
input_mod1 = ad.read_h5ad(testpar['input_mod1'])
output = ad.read_h5ad(testpar['output'])

print("> Checking contents of output.h5ad")
assert output.uns['dataset_id'] == input_mod1.uns['dataset_id']
assert output.uns['method_id'] == method_id
assert output.n_obs == input_mod1.n_obs
assert output.n_vars >= 1
assert output.n_vars <= 100
assert all(output.obs_names == input_mod1.obs_names)

print("> Test succeeded!")