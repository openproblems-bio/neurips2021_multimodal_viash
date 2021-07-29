import anndata
import scipy.spatial
from scipy.sparse import csr_matrix
import numpy as np

# VIASH START
par = {
    "input_mod1": "../../../../resources_test/predict_modality/test_resource.mod1.h5ad",
    "input_mod2": "../../../../resources_test/predict_modality/test_resource.mod2.h5ad",
    "output": "../../../../resources_test/predict_modality/test_resource.prediction.h5ad",
}
# VIASH END

# load dataset to be censored
ad_rna = anndata.read_h5ad(par["input_mod1"])
ad_mod2 = anndata.read_h5ad(par["input_mod2"])


# Find the correct shape

ad_train = ad_rna[ad_rna.obs["group"] == "train"]
ad_test = ad_rna[ad_rna.obs["group"] == "test"]

# Note this is to test spare prediction matricies, should this be allowed?
prediction = csr_matrix((ad_test.X.shape[0], ad_mod2.shape[1]), dtype = np.float32)

# Write out prediction
out = anndata.AnnData(
    X=prediction,
    uns={
        "dataset_id": ad_rna.uns["dataset_id"],
        "method_id": "baseline_zeros",
    }
)
out.write_h5ad(par["output"])
