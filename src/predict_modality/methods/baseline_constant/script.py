import anndata
import scipy.spatial
from scipy.sparse import csr_matrix
import numpy as np

# VIASH START
par = {
    "input_mod1_train": "../../../../resources_test/predict_modality/test_resource.mod1.h5ad",
    "input_mod1_test": "../../../../resources_test/predict_modality/test_resource.mod1.h5ad",
    "input_mod2": "../../../../resources_test/predict_modality/test_resource.mod2.h5ad",
    "output": "../../../../resources_test/predict_modality/test_resource.prediction.h5ad",
}
# VIASH END

# load dataset to be censored
ad_rna_train = anndata.read_h5ad(par["input_mod1_train"])
ad_rna_test = anndata.read_h5ad(par["input_mod1_test"])
ad_mod2 = anndata.read_h5ad(par["input_mod2"])


# Find the correct shape

ad_mod2_train = ad_mod2[ad_mod2.obs["group"] == "train"]
mean = np.array(ad_mod2_train.X.mean(axis=0)).flatten()
prediction = np.tile(mean, (ad_rna_test.shape[0], 1))

# Write out prediction
out = anndata.AnnData(
    X=prediction,
    uns={
        "dataset_id": ad_rna.uns["dataset_id"],
        "method_id": "baseline_zeros",
    }
)
out.write_h5ad(par["output"])
