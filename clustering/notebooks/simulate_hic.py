import numpy as np

def tads(k, distance, alpha):
    interaction = k*distance**alpha
    return interaction

def _sim_mat(n, k, alpha):
    cell1 = np.zeros((n, n))
    
    ij = np.ones((n, n))*np.linspace(0, n-1, n)
    i = ij.T[np.triu_indices(n)]
    j = ij[np.triu_indices(n)]
    
    distance = j - i + 1
    
    Bmean = tads(k, distance, alpha)
    
    cell1[np.triu_indices(n)] = Bmean
    cell1 = cell1 + cell1.T
    cell1[np.diag_indices(n)] /= 2
    
    return cell1

class SimulateHiC: