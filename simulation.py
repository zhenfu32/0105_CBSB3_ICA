import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
from collections import Counter

def oligo_pool(pool_size,k):
    """
    This function can create the origin oligonucleotides pool.
    :param n: pool size
    :param k: the length of sequence
    :return: the origin pool
    """
    sequences = []
    bases = "ATCG"
    for i in range(pool_size):
        seq = "".join(random.choices(bases, k=k))
        sequences.append(seq)
    return sequences

def simulate_selex(pool_size, rounds, affinity):
    """
    the main simulation function
    :param pool_size: the pool size of each rounds
    :param rounds:
    :param affinity: the binding affinity matrix of target
    :return:
    """
    pool = np.zeros(shape=(pool_size, rounds))
    for i in range(rounds):
       # normalize the affinity values
       norm_affinity = affinity / np.sum(affinity)
       # select oligonucleotides based on affinity
       selected_indices = np.random.choice(pool_size, size=pool_size, replace=True, p= norm_affinity)
       selected_pool = np.zeros(shape=(pool_size, rounds))
       selected_pool[selected_indices, i] = 1
       pool[:, i] = np.sum(selected_pool, axis=1)
       # update the affinity based on the selected oligonucleotides
       affinity[selected_indices] = np.random.normal(loc=affinity[selected_indices], scale=0.1)
    return pool



# The first step: create oligonucleotides pool
# In this simulation, we create the pool with 1 million random 16mer sequences
pool_size = 1000000
k_simulation = 16
rounds = 3
pool_r0 = oligo_pool(pool_size, k_simulation)
# print(pool_r0[:5])
# count the amount of sequences
count_r0 = Counter(pool_r0)
# print(count_r0.most_common(10),len(count_r0))


# The next step: import Gibbs energy mateix estimated from a SELEX experiment
# on the transcription factor Bicoid
gibbs_matrix = np.array([
    [-4.722516, -5.729347, 0.000000, -6.251779],
    [-7.447426, -5.981440, 0.000000, -16.853690],
    [0.000000, -6.946246, -15.701235, -8.529272],
    [-7.746046, -15.548042, -12.535315, 0.000000],
    [-7.989755, -7.201358, -24.708969, 0.000000],
    [0.000000, -9.611195, -8.497223, -5.336888],
    [-0.505663, -19.926999, 0.000000, -4.445374],
    [-1.836787, -0.228140, 0.000000, -0.945140],
    [-1.841359, -1.612913, 0.000000, -1.417988],
    [-1.431632, -1.539663, 0.000000, -0.235633]])


affinity = np.random.normal(loc=1, scale=0.1, size=pool_size)
pool = simulate_selex(pool_size, rounds, affinity)
print(pool[:5,])

# create a DataFrame to store the SELEX data
df = pd.DataFrame(data=pool, columns=[f'Round {i + 1}' for i in range(rounds)])
print(df[:5])
# # compute the biochemical parameters based on the SELEX data
parameters = df.apply(lambda x: np.log(x) - np.log(pool_size - x)).mean()
# # compute the affinity of each oligonucleotide
affinity_estimate = np.exp(np.sum(np.outer(affinity, parameters), axis=1))
