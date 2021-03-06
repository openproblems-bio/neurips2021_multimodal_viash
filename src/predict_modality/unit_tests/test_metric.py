from os import path
import subprocess
import anndata as ad
import pandas as pd

## VIASH START
# This code block will be replaced by viash at runtime.
meta = { 'functionality_name': 'foo' }
meta_path = "src/predict_modality/metrics/check_format/metric_meta_check_format.tsv"
## VIASH END

method_id = meta['functionality_name']
command = "./" + method_id

# define some filenames
testpar = {
  "input_prediction": "resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
  "input_solution": "resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod2.h5ad",
  "output": "output.h5ad"
}
meta_path = resources_dir + '/metric_meta.tsv'

print("> Running method")
out = subprocess.check_output([
  command,
  "--input_prediction", testpar['input_prediction'],
  "--input_solution", testpar['input_solution'],
  "--output", testpar['output']
]).decode("utf-8")

print("> Checking whether output files were created")
assert path.exists(testpar['output'])

print("> Reading h5ad files")
input_prediction = ad.read_h5ad(testpar['input_prediction'])
input_solution = ad.read_h5ad(testpar['input_solution'])
output = ad.read_h5ad(testpar['output'])

metric_meta = pd.read_csv(
  meta_path, 
  delimiter="\t",
  header=0,
  dtype={ 'metric_id': str, 'metric_min': float, 'metric_max': float, 'metric_higherisbetter': bool }
)

print("> Checking contents of metric_meta.tsv")
assert 'metric_id' in metric_meta
assert 'metric_min' in metric_meta
assert 'metric_max' in metric_meta
assert 'metric_higherisbetter' in metric_meta

print("> Checking .uns['dataset_id']")
assert 'dataset_id' in output.uns
assert output.uns['dataset_id'] == input_prediction.uns['dataset_id']

print("> Checking .uns['method_id']")
assert 'method_id' in output.uns
assert output.uns['method_id'] == input_prediction.uns['method_id']

print("> Checking .uns['metric_ids']")
assert 'metric_ids' in output.uns
assert set(output.uns['metric_ids']) == set(metric_meta.metric_id)

print("> Checking .uns['metric_values']")
assert 'metric_values' in output.uns
assert output.uns['metric_ids'].size == output.uns['metric_values'].size

# merge with metric_meta to see if metric_value lies within the expected range
output_uns = pd.DataFrame({
  'metric_id': output.uns['metric_ids'], 
  'metric_value': output.uns['metric_values']
})

scores = metric_meta.merge(output_uns, on="metric_id")

assert all(scores.metric_value >= scores.metric_min)
assert all(scores.metric_value <= scores.metric_max)

print("> Test succeeded!")