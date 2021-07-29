## VIASH START
par = dict(
    input_prediction="resources_test/joint_embedding/test_resource.prediction.h5ad",
    input_solution="resources_test/joint_embedding/test_resource.solution.h5ad",
    output="resources_test/joint_embedding/test_resource.asw_batch.tsv",
    debug=True
)

## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
from scIB.metrics import silhouette_batch

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'asw_batch'

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

print('Preprocessing')
adata.obsm['X_emb'] = adata.X
sc.pp.neighbors(adata, use_rep='X_emb')

print('Compute score')
_, sil_clus = silhouette_batch(
    adata,
    batch_key='batch',
    group_key='cell_type',
    embed='X_emb',
    verbose=False
)
score = sil_clus['silhouette_score'].mean()

with open(output, 'w') as file:
    header = ['dataset', 'output_type', 'metric', 'value']
    entry = [dataset_id, OUTPUT_TYPE, METRIC, score]
    file.write('\t'.join(header) + '\n')
    file.write('\t'.join([str(x) for x in entry]))
