import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
from collections import Counter

def oligo_pool(n,k):
    sequences = []
    bases = "ATCG"
    for i in range(n):
        seq = "".join(random.choices(bases, k=k))
        sequences.append(seq)
    return sequences




# The first step: create oligonucleotides pool
# In this simulation, we create the pool with 1 million random 16mer sequences
n_simulation = 1000000
k_simulation = 16
pool_r0 = oligo_pool(n_simulation, k_simulation)
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
