import cooler

import numpy as np


class DataPreparation:
    def __init__(self, filename, resolution, region, balance=True):
        """
        Transform HiC data on .mcool format to 2D numpy arrays for further analysis.

        Args:
            filename (string):  
                Path to a multi-resolution coolers .mcool file
            resolution (int): 
                Resolution for cooler object
            region  (string):
                (Genomic range strings may be provided as 3-tuples (chrom: str, start: int, end: int) or in UCSC-style strings of the style {chrom}:{start}-{end}. Unit prefixes k, M, G are supported in range strings. For regions with start and end that are not multiples of the resolution, selectors return the range of shortest range bins that fully contains the open interval [start, end)) 
                Genomic range string of the style {chrom}:{start}-{end}, unit 
                prefixes k, M, G are supported
            balance (bool, optional):
                default=True, whether to balance data or not
        Attr:
            matrix (ndarray):
                Matrix of HiC data for given resolution and region
        """
        self.filename = filename
        self.resolution = resolution
        self.region = region
        self.clr, self.matrix = self.create_matrix(balance)
        
    def create_matrix(self, balance=True):
        """
        Loads Hi-C map at a given resolution from a cooler file

        Args:
            balance (bool, optional):
                default=True, whether to balance data or not

        Return:
            matrix (ndarray):
                Matrix of HiC data for given resolution and region
        """
        clr = cooler.Cooler(self.filename
                            + '::resolutions/'
                            + str(self.resolution))
        matrix = clr.matrix(balance=balance).fetch(self.region)
        return clr, matrix
    
    def higlass_ratio(self, other, replace_zero_zero=False):
        """Calculates the difference between two Hi-C maps by division.
        
        The calculations are done in the same way as in HiGlass' divide by.
        All division by zero is set to NaN.
        
        Args:
            other (class instance): 
                Instance of class DataPreparation
            replace_zero_zero (bool, optional):
                default=False, Replace all 0/0 by 1

        Return:
            diff (array):
                The two matrices divided
        """
        if not hasattr(self, 'matrix'):
            self.create_matrix()
        elif not hasattr(other, 'matrix'):
            other.create_matrix()
            
        # replaces all zeros in denominator to NaN to ensure no division by zero
        other.matrix[other.matrix == 0] = np.nan
        
        diff = self.matrix/other.matrix
        
        # add pseudocount - the smallest value in diff
        #diff += np.nanmin(diff)
        
        if replace_zero_zero:
            # find all 0/0
            sum_data = self.matrix + other.matrix

            # set all 0/0 to 1
            diff[sum_data == 0] = 1
        
        return diff
    
    def divide_with_pseudocount(self, other, pseudocount=1):
        pseudo_mat = np.ones(self.matrix.shape)*pseudocount
        
        diff = (self.matrix + pseudo_mat) / (other.matrix + pseudo_mat)
        
        return diff
    
    def subtract(self, other):
        """
        Calculates the difference between two Hi-C maps by subtraction.

        Args:
            other (class instance): 
                Instance of class DataPreparation

        Return:
            diff (array):
                The two matrices divided
        """
        if not hasattr(self, 'matrix'):
            self.create_matrix()
        elif not hasattr(other, 'matrix'):
            other.create_matrix()
            
        diff = self.matrix - other.matrix
        return diff
    
    def relative_difference(self, other):
        """Calculates the relative difference between two Hi-C contact matrices.
        
        Let x and y be pairwise matrix elements. Then the relative difference
        is calculated by
            (x - y)/( (x + y)/2 ).
        x and y are equal when the relative difference is 0. When both x and y
        are 0, the mean ( (x + y)/2 ) is set to 1 to ensure the result is 0.
        
        Args:
            other (class instance): 
                Instance of class DataPreparation

        Return:
            diff (array):
                The two matrices divided
        """
        if not hasattr(self, 'matrix'):
            self.create_matrix()
        elif not hasattr(other, 'matrix'):
            other.create_matrix()
        
        # calculate the average IF
        mean = (self.matrix + other.matrix)/2
        
        # make sure that pairs where both IF are zero returns 0
        mean[mean == 0] = 1
        
        # calculate the relative difference
        diff = (self.matrix - other.matrix)/mean
        return diff

    def raw_data(self):
        """
        Creates numpy array containing raw data for given region.

        Return:
            raw_data (array): numpy array containing the raw data for a given
                              region
        """
        raw_data = self.clr.matrix(balance=False).fetch(self.region)
        return raw_data

    def balanced_data(self):
        """
        Creates numpy array containing the balanced data for a given region.

        Return:
            balanced_data (array): numpy array containing the balanced data for
                                   a given region
        """
        balanced_data = self.clr.matrix(balance=True).fetch(self.region)
        return balanced_data



if __name__ == '__main__':
    filename = '/home/elinfi/coolers/HiC_wt_001.mcool'
    resolution = 1024000
    region = 'chr4'
    
    test = DataPreparation(filename, resolution, region)
    test2 = DataPreparation(filename, resolution, region)
    test.create_matrix()
    test.divide(test2)
