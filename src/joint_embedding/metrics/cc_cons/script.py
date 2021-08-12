## VIASH START
par = dict(
    input_prediction="resources_test/joint_embedding/test_resource.prediction.h5ad",
    input_solution="resources_test/joint_embedding/test_resource.solution.h5ad",
    output="resources_test/joint_embedding/test_resource.cc_cons.tsv",
    debug=True
)
## VIASH END

debug = par['debug']

print('Importing libraries')
import pprint
import scanpy as sc
import anndata
import numpy as np
from scIB.metrics import cell_cycle

if debug:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'cc_cons'

input_prediction = par['input_prediction']
input_solution = par['input_solution']
output = par['output']

print("Read prediction anndata")
adata = sc.read(input_prediction)
dataset_id = adata.uns['dataset_id']

if 'organism' in adata.uns:
    organism = adata.uns['organism']
else:
    organism = 'synthetic'

print("Read solution anndata")
adata_solution = sc.read(input_solution)
if debug:
    print(adata_solution.X)
    print(adata_solution.var)
    print(adata.var)

print('Transfer obs annotations')
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs_names]
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs_names]
recompute_cc = not np.all(np.in1d(['S_score', 'G2M_score'], adata_solution.obs_keys()))

print('Preprocessing')
adata.obsm['X_emb'] = adata.X

print('Compute score')
if organism == 'synthetic':
    score = np.nan if recompute_cc else cell_cycle(
        adata_pre=adata_solution,
        adata_post=adata,
        batch_key='batch',
        embed='X_emb',
        recompute_cc=False
    )
else:
    score = cell_cycle(
        adata_pre=adata_solution,
        adata_post=adata,
        batch_key='batch',
        embed='X_emb',
        recompute_cc=recompute_cc,
        organism=organism
    )

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

# # store score as tsv
# with open(output, 'w') as file:
#     header = ['dataset', 'output_type', 'metric', 'value']
#     entry = [dataset_id, OUTPUT_TYPE, METRIC, score]
#     file.write('\t'.join(header) + '\n')
#     file.write('\t'.join([str(x) for x in entry]))
