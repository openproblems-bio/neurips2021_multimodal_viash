
import numpy as np
import scipy
from sklearn.neighbors import NearestNeighbors

def entropy_batch_mixing(
    latent_space, batches, n_neighbors=50, n_pools=50, n_samples_per_pool=100
):
    def entropy(hist_data):
        n_batches = len(np.unique(hist_data))
        if n_batches > 2:
            raise ValueError("Should be only two clusters for this metric")
        frequency = np.mean(hist_data == 1)
        if frequency == 0 or frequency == 1:
            return 0
        return -frequency * np.log(frequency) - (1 - frequency) * np.log(1 - frequency)

    def neg_kl(hist_data, global_freq):
        n_batches = len(np.unique(hist_data))
        if n_batches > 2:
            raise ValueError("Should be only two clusters for this metric")
        frequency = np.mean(hist_data == 1)
        if frequency == 0 or frequency == 1:
            return 0
        return -(
            frequency * np.log(frequency / global_freq)
            + (1 - frequency) * np.log((1 - frequency) / (1 - global_freq))
        )

    n_neighbors = min(n_neighbors, len(latent_space) - 1)
    nne = NearestNeighbors(n_neighbors=1 + n_neighbors, n_jobs=8)
    nne.fit(latent_space)
    kmatrix = nne.kneighbors_graph(latent_space) - scipy.sparse.identity(
        latent_space.shape[0]
    )

    global_freq = np.mean(batches)
    print(global_freq)
    score = 0
    for t in range(n_pools):
        indices = np.random.choice(
            np.arange(latent_space.shape[0]), size=n_samples_per_pool
        )
        score += np.mean(
            [
                neg_kl(
                    batches[  # the batches of cell i's neighbors
                        kmatrix[indices].nonzero()[
                            1
                        ][  # the neighbors of cell i (columns in row i)
                            kmatrix[indices].nonzero()[0] == i  # the row of cell i
                        ]
                    ],
                    global_freq,
                )
                for i in range(n_samples_per_pool)
            ]
        )
    return score / float(n_pools)
