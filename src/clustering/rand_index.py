import numpy as np

from sklearn.metrics.cluster import rand_score, adjusted_rand_score

def rand_index(list_filenames):
    """
    Args:
        list_filenames (list of string):
            List of files containing cluster labels to be compared.
    """
    initial_data = np.load(list_filenames[0])
    
    n = initial_data.shape[0]
    triu_idx = np.triu_indices(n)
    nan_idx = np.isnan(initial_data)
    nan_idx_triu = nan_idx[triu_idx]
    print(nan_idx_triu.shape)
    
    m = len(list_filenames)
    rand_idxs = np.zeros((m, m))
    
    for i in range(m):
        for j in range(i, m):
            print(np.unique(np.load(list_filenames[i])))
            mat1 = np.load(list_filenames[i])[triu_idx]#.astype(int)
            mat2 = np.load(list_filenames[j])[triu_idx]#.astype(int)
            
            # remove all nan's in both
            nan1 = np.isnan(mat1)
            nan2 = np.isnan(mat2)
            nans = nan1 + nan2
            mat1 = mat1[~nans]
            mat2 = mat2[~nans]
            
            print(np.unique(mat1))
            print(np.unique(mat2))
            
            rand_idx = adjusted_rand_score(mat1, mat2)
            rand_idxs[i, j] = rand_idx
            rand_idxs[j, i] = rand_idx
    
    return rand_idxs
