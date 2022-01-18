import numpy as np

from data_preparation import DataPreparation

def _powerlaw(distance, C, alpha):
    result = np.where(distance == 0, C, C * distance**(-alpha))
    return result

def _poisson(distance):
    result = np.zeros(distance)
    
    for i, lmbda in enumerate(distance):
        result[i] = np.random.Generator.poisson(lmbda)
    return result

def _sim_mat(n, medianIF, powerlaw_alpha):
    cell1 = np.zeros((n, n))
    cell2 = np.zeros((n, n))
    
    ij = np.ones((n, n))*np.linspace(0, n-1, n)
    i = ij.T[np.triu_indices(n)]
    j = ij[np.triu_indices(n)]
    
    distance = j - i + 1
    
    #Bmean = _powerlaw(distance, medianIF, powerlaw_alpha)
    Bmean = _poisson(distance)
    
    cell1[np.triu_indices(n)] = Bmean
    cell1 = cell1 + cell1.T
    cell1[np.diag_indices(n)] /= 2
    
    return cell1


class SimulateData():
    def __init__(self, resolution, region):
        """
        Simulate Hi-C data.
        
        Args:
            resolution (int):
                Resolution for cooler object
            region  (string): 
                Genomic range string of the style {chrom}:{start}-{end}, unit 
                prefixes k, M, G are supported
        Attrs:
            matrix (ndarray):
                Numpy zeroes array.
        """
        self.resolution = resolution
        self.region = region
        self.clr, self.matrix = self.initial_matrix(resolution, region)
        
    def initial_matrix(self, resolution, region):
        """
        Creates a matrix with background noise.
        
        Args:
            resolution (int):
                Resolution for cooler object
            region  (string): 
                Genomic range string of the style {chrom}:{start}-{end}, unit 
                prefixes k, M, G are supported
                
        Returns:
            matrix (ndarray):
                Matrix with the remaining noise difference from the two 
                replicates wild type 001 and wild type 002.
        """
        path_wt_001 = '/home/elinfi/coolers/HiC_wt_001.mcool'
        path_wt_002 = '/home/elinfi/coolers/HiC_wt_002.mcool'
        
        wt_001 = DataPreparation(path_wt_001, resolution, region)
        wt_002 = DataPreparation(path_wt_002, resolution, region)
        
        matrix = wt_001.matrix - wt_002.subtract(wt_001)
        
        return wt_001.clr, matrix
    
    def addTAD(self, tad_mat, ul):
        test = np.zeros_like(self.matrix)
        test[ul:ul + tad_mat.shape[0], ul:ul + tad_mat.shape[1]] = tad_mat
        
        return self.matrix + test
    
    
        
    
    
        
    
        