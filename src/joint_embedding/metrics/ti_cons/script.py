## VIASH START
par = dict(
    input_prediction="resources_test/joint_embedding/test_resource.prediction.h5ad",
    input_solution="resources_test/joint_embedding/test_resource.solution.h5ad",
    output="resources_test/joint_embedding/test_resource.ti_cons.tsv",
    debug=True
)
## VIASH END

debug = par['debug']

print('Importing libraries')
import pprint
import numpy as np
import scanpy as sc
import anndata
from scIB.metrics import trajectory_conservation

if debug:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'ti_cons'

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = sc.read(input_prediction)
dataset_id = adata.uns['dataset_id']

print("Read solution anndata")
adata_solution = sc.read(input_solution)

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs_names]
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs_names]
adt_atac_trajectory = 'ATAC_trajectory'

print('Compute scores')
obs_keys = adata_solution.obs_keys()

if 'RNA_trajectory' in obs_keys:
    score_rna = trajectory_conservation(
        adata_pre=adata_solution,
        adata_post=adata,
        label_key='cell_type',
        pseudotime_key='RNA_trajectory'
    )
else:
    score_rna = np.nan

if 'ADT_trajectory' in obs_keys or 'ATAC_trajectory' in obs_keys:
    score_adt_atac = trajectory_conservation(
        adata_pre=adata_solution,
        adata_post=adata,
        label_key='cell_type',
        pseudotime_key=adt_atac_trajectory
    )
else:
    score_adt_atac = np.nan

score_mean = (score_rna + score_adt_atac) / 2

# store adata with metrics
print("Create output object")
out = anndata.AnnData(
    uns=dict(
        dataset_id=adata.uns['dataset_id'],
        method_id=adata.uns['method_id'],
        metric_ids=['ti_cons_RNA', 'ti_cons_ADT_ATAC', 'ti_cons_mean'],
        metric_values=[score_rna, score_adt_atac, score_mean]
    )
)

print("Write output to h5ad file")
out.write(output, compression='gzip')

# # store score as tsv
# with open(output, 'w') as file:
#     header = ['dataset', 'output_type', 'metric', 'value']
#     entry = [dataset_id, OUTPUT_TYPE, METRIC, score]
#     file.write('\t'.join(header) + '\n')
#     file.write('\t'.join([str(x) for x in entry]))
