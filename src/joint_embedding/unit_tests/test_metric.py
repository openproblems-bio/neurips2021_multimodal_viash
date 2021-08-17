from os import path
import subprocess
import anndata as ad
import pandas as pd

par = {
  "input_prediction": "resources_test/joint_embedding/test_resource.prediction.h5ad",
  "input_solution": "resources_test/joint_embedding/test_resource.solution.h5ad",
  "output": "output.h5ad",
  "metric_meta": "metric_meta.tsv"
}

print("> Running method")
out = subprocess.check_output([
  "./" + meta['functionality_name'],
  "--input_prediction", par['input_prediction'],
  "--input_solution", par['input_solution'],
  "--output", par['output']
]).decode("utf-8")

print("> Checking whether output files were created")
assert path.exists(par['output'])

print("> Reading h5ad files")
input_prediction = ad.read_h5ad(par['input_prediction'])
input_solution = ad.read_h5ad(par['input_solution'])
output = ad.read_h5ad(par['output'])

metric_meta = pd.read_csv(
  resources_dir + '/metric_meta.tsv',
  delimiter="\t"
)

print("> Checking contents of metric_meta.tsv")
assert 'metric_id' in metric_meta
assert 'metric_min' in metric_meta
assert 'metric_max' in metric_meta
assert 'metric_higherisbetter' in metric_meta

print("> Checking contents of output.h5ad")
assert 'dataset_id' in output.uns
assert output.uns['dataset_id'] == input_prediction.uns['dataset_id']

assert 'method_id' in output.uns
assert output.uns['method_id'] == input_prediction.uns['method_id']

assert 'metric_ids' in output.uns
assert 'metric_values' in output.uns
assert set(output.uns['metric_ids']) == set(metric_meta.metric_id)

score = pd.DataFrame({
  'metric_id': output.uns['metric_ids'], 
  'metric_value': output.uns['metric_values']
})

comb = metric_meta.merge(score, on="metric_id")

assert all(comb.metric_value >= comb.metric_min)
assert all(comb.metric_value <= comb.metric_max)

print("> Test succeeded!")