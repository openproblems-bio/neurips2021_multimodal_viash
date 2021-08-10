print("Loading dependencies")
import anndata as ad
import random
import numpy as np
import scipy.sparse

# VIASH START
par = {
    "input_mod1": "resources_test/common/test_resource.output_rna.h5ad",
    "input_mod2": "resources_test/common/test_resource.output_mod2.h5ad",
    "output_train_mod1": "resources_test/match_modality/test_resource.train_mod1.h5ad",
    "output_train_mod2": "resources_test/match_modality/test_resource.train_mod2.h5ad",
    "output_train_sol": "resources_test/match_modality/test_resource.train_sol.h5ad",
    "output_test_mod1": "resources_test/match_modality/test_resource.test_mod1.h5ad",
    "output_test_mod2": "resources_test/match_modality/test_resource.test_mod2.h5ad",
    "output_test_sol": "resources_test/match_modality/test_resource.test_sol.h5ad"
}
# VIASH END

print("Reading input data")
mod1 = ad.read_h5ad(par["input_mod1"])
mod2 = ad.read_h5ad(par["input_mod2"])
new_dataset_id = mod1.uns["dataset_id"] + "_MM"

print("Shuffling train cells")
train_mod1_ix = np.where(mod1.obs["is_train"])[0]
train_mod2_ix = train_mod1_ix.copy()
random.shuffle(train_mod2_ix)
train_mod1 = mod1[train_mod1_ix, :]
train_mod2 = mod2[train_mod2_ix, :]

print("Shuffling test cells")
test_mod1_ix = np.where(~mod1.obs["is_train"])[0]
test_mod2_ix = test_mod1_ix.copy()
random.shuffle(test_mod2_ix)
test_mod1 = mod1[test_mod1_ix, :]
test_mod2 = mod2[test_mod2_ix, :]


print("Creating mod1 outputs")
desired_var1_cols = [x for x in ["gene_ids", "feature_types"] if x in mod1.var.columns]
out_train_mod1 = ad.AnnData(
    X=train_mod1.X,
    var=mod1.var[desired_var1_cols],
    uns={ "dataset_id": new_dataset_id },
)
out_train_mod1.X.sort_indices()
out_test_mod1 = ad.AnnData(
    X=test_mod1.X,
    var=mod1.var[desired_var1_cols],
    uns={ "dataset_id": new_dataset_id },
)
out_test_mod1.X.sort_indices()

print("Creating mod2 outputs")
desired_var2_cols = [x for x in ["gene_ids", "feature_types"] if x in mod2.var.columns]
out_train_mod2 = ad.AnnData(
    X=train_mod2.X,
    var=mod2.var[desired_var2_cols],
    uns={ "dataset_id": new_dataset_id },
)
out_train_mod2.X.sort_indices()
out_test_mod2 = ad.AnnData(
    X=test_mod2.X,
    var=mod2.var[desired_var2_cols],
    uns={ "dataset_id": new_dataset_id },
)
out_test_mod2.X.sort_indices()

print("Creating solution outputs")
out_train_sol_mat = scipy.sparse.csr_matrix(
  (np.ones(train_mod1_ix.size), (train_mod1_ix, train_mod2_ix))
)
out_train_sol = ad.AnnData(
    X=out_train_sol_mat,
    uns={ "dataset_id": new_dataset_id },
    dtype="float32",
)
out_test_sol_mat = scipy.sparse.csr_matrix(
  (np.ones(test_mod1_ix.size), (test_mod1_ix, test_mod2_ix))
)
out_test_sol = ad.AnnData(
    X=out_test_sol_mat,
    uns={ "dataset_id": new_dataset_id },
    dtype="float32"
)

print("Writing output objects to file")
out_train_mod1.write_h5ad(filename=par["output_train_mod1"], compression="gzip")
out_train_mod2.write_h5ad(filename=par["output_train_mod2"], compression="gzip")
out_train_sol.write_h5ad(filename=par["output_train_sol"], compression="gzip")
out_test_mod1.write_h5ad(filename=par["output_test_mod1"], compression="gzip")
out_test_mod2.write_h5ad(filename=par["output_test_mod2"], compression="gzip")
out_test_sol.write_h5ad(filename=par["output_test_sol"], compression="gzip")
