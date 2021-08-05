from os import path
import subprocess
import anndata as ad
# import pandas as pd
import numpy as np

np.random.seed(42)

metric = 'ari'
# metric_file = metric + '.tsv'
metric_file = metric + '.h5ad'

print(">> Running script")
out = subprocess.check_output([
    "./" + metric,
    "--input_prediction", 'resources_test/joint_embedding/test_resource.prediction.h5ad',
    "--input_solution", 'resources_test/joint_embedding/test_resource.solution.h5ad',
    "--output", metric_file
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(metric_file)
# result = pd.read_table(metric_file)
result = ad.read_h5ad(metric_file)
sol = ad.read_h5ad('resources_test/joint_embedding/test_resource.solution.h5ad')
pred = ad.read_h5ad('resources_test/joint_embedding/test_resource.prediction.h5ad')

# print(">> Check that score makes sense")
# assert result.shape == (1, 4)
# score = result.loc[0, 'value']

print(">> Check contents of result.uns")
assert 'dataset_id' in result.uns
assert result.uns['dataset_id'] == sol.uns['dataset_id']

assert 'method_id' in result.uns
assert result.uns['method_id'] == pred.uns['method_id']

assert 'metric_ids' in result.uns
assert result.uns['metric_ids'] == [metric]

assert 'metric_values' in result.uns
score = result.uns['metric_values'][0]
print(score)
assert 0 <= score <= 1
assert score == 0.336331837241985

print(">> All tests passed successfully")
