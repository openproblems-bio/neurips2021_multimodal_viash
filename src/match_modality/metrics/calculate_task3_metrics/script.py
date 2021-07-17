import math

import anndata
import numpy as np

# VIASH START
par = {
    "prediction": "../../../../resources_test/common/pbmc_1k_protein_v3.prediction.h5ad",
    "solution": "../../../../resources_test/common/pbmc_1k_protein_v3.solution.h5ad",
}
# VIASH END

# load predictions and solutions
ad_prediction = anndata.read_h5ad(par["prediction"])
ad_solution = anndata.read_h5ad(par["solution"])

prediction = ad_prediction.uns["prediction"]
solution = ad_solution.uns["pairings"].toarray()

prediction = prediction.flatten(order="C")
solution = solution.flatten(order="C")

ordering = np.argsort(prediction)[::-1]

sorted_prediction = np.around(prediction[ordering])
sorted_solution = solution[ordering]
true_predictions = sorted_prediction == sorted_solution

# Calculate metrics
num_selected = np.array([i for i in range(len(true_predictions))])
tp = np.cumsum(true_predictions)
fp = num_selected - tp
num_possible_pairings = len(solution)
num_positive_pairings = sum(solution)
num_negative_pairings = len(solution) - sum(solution)

fn = num_positive_pairings - tp
tn = num_negative_pairings - fp
acc = (tp + tn) / (num_positive_pairings +num_negative_pairings)
tpr = tp / num_positive_pairings
spec = tn / num_negative_pairings
prec = 1 if len(num_selected) == 0 else tp / (tp + fp)
npv = tn / (tn + fn)
f1 = 2 * tp / (2 * tp + fp + fn)
# mcc = 0 if len(num_selected) == 0 else (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
informedness = tpr + spec - 1
markedness = prec + npv - 1

auroc = np.trapz(1 - spec, tpr)
aupr = abs(np.trapz(tpr, prec))
f1 = 2 * auroc * aupr / (auroc + aupr) if auroc + aupr != 0 else 0

thing = 0
