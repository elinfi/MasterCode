import numpy as np

from data_preparation import DataPreparation

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
        
        matrix = wt_002.subtract(wt_001)
        
        return wt_001.clr, matrix
    
    
        
    
        