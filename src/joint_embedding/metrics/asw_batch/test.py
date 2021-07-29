from os import path
import subprocess
import pandas as pd
import numpy as np

np.random.seed(42)

metric = 'asw_batch'
metric_file = metric + '.tsv'

print(">> Running script")
out = subprocess.check_output([
    "./" + metric,
    "--input_prediction", 'resources_test/joint_embedding/test_resource.prediction.h5ad',
    "--input_solution", 'resources_test/joint_embedding/test_resource.solution.h5ad',
    "--output", metric_file
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(metric_file)
result = pd.read_table(metric_file)

print(">> Check that score makes sense")
assert result.shape == (1, 4)
assert result['metric'][0] == metric
score = result.loc[0, 'value']
print(score)

assert 0 < score < 1
assert score == 0.7223580315216693

print(">> All tests passed successfully")
