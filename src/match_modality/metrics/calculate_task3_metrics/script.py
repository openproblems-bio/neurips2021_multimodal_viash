import anndata
import scipy.sparse
import numpy as np

from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

# VIASH START
par = {
    "prediction": "../../../../prediction.h5ad",
    "solution": "../../../../prediction.h5ad"
}
# VIASH END

# load dataset to be censored
prediction = anndata.read_h5ad(par["prediction"])
solution = anndata.read_h5ad(par["solution"])

# Get predicted pairings & rankings
prediction_matrix = prediction.obs["prediction"]

# Get actual pairings
true_pairing = solution.uns["solution"]




