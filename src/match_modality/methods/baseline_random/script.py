import anndata
import scipy.spatial
import numpy as np

# VIASH START
par = {
    "input_mod1": "../../../../resources_test/match_modality/test_resource.mod1.h5ad",
    "input_mod2": "../../../../resources_test/match_modality/test_resource.mod2.h5ad",
    "output": "../../../../resources_test/match_modality/test_resource.prediction.h5ad",
}
# VIASH END

# load dataset to be censored
ad_rna = anndata.read_h5ad(par["input_mod1"])
ad_mod2 = anndata.read_h5ad(par["input_mod2"])

rna, mod2 = ad_rna.X, ad_mod2.X

pairing_matrix = np.zeros((rna.shape[0], mod2.shape[0]))
indices = np.random.randint(mod2.shape[0], size=rna.shape[0])
pairing_matrix[np.arange(rna.shape[0]), indices] = 1

# Write out prediction
prediction = anndata.AnnData(
    X=pairing_matrix,
    uns={
        "dataset_id": ad_rna.uns["dataset_id"],
    }
)
prediction.write_h5ad(par["output"])
