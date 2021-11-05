## VIASH START
par = dict(
    input_prediction="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
    input_solution="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.solution.h5ad",
    output="resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.graph_conn.tsv",
    debug=True
)

## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
import anndata
import scib

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'graph_conn'

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
score = scib.me.graph_connectivity(adata, label_key='cell_type')

# store adata with metrics
print("Create output object")
out = anndata.AnnData(
    uns=dict(
        dataset_id=adata.uns['dataset_id'],
        method_id=adata.uns['method_id'],
        metric_ids=[METRIC],
        metric_values=[score]
    )
)

print("Write output to h5ad file")
out.write(output, compression='gzip')

# store score as tsv
# with open(output, 'w') as file:
#     header = ['dataset', 'output_type', 'metric', 'value']
#     entry = [dataset_id, OUTPUT_TYPE, METRIC, score]
#     file.write('\t'.join(header) + '\n')
#     file.write('\t'.join([str(x) for x in entry]))
