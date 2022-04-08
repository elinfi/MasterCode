import numpy as np

from sklearn.metrics.cluster import rand_score

def rand_index(list_filenames):
    """
    Args:
        list_filenames (list of string):
            List of files containing cluster labels to be compared.
    """
    n = len(list_filenames)
    rand_idxs = np.zeros((n, n))
    triu_idx = np.triu_indices(n)
    
    for i in range(n):
        for j in range(i, n):
            mat1 = np.load(list_filenames[i])[triu_idx]
            mat2 = np.load(list_filenames[j])[triu_idx]
            rand_idx = rand_score(mat1, mat2)
            rand_idxs[i, j] = rand_idx
            rand_idxs[j, i] = rand_idx
    
    return rand_idxs