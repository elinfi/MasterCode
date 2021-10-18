import numpy as np

class SimulateData:
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
        self.resolutioin = resolution
        self.region = region
        self.matrix = self.create_matrix(resolution, region)
        
    def create_matrix(self, resolution, region):
        """
        Creates empty matrix filled with zeroes.
        
        Args:
            resolution (int):
                Resolution for cooler object
            region  (string): 
                Genomic range string of the style {chrom}:{start}-{end}, unit 
                prefixes k, M, G are supported
                
        Returns:
            matrix (ndarray):
                Zeroes numpy array.
        """
        n = 50
        matrix = np.zeros((n, n))
        