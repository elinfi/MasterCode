"""Methods for comparing two Hi-C matrices"""

import numpy as np

def relative_difference(mat1, mat2):
    """Calculates the relative difference between two Hi-C contact matrices.
        
        Let x and y be pairwise matrix elements. Then the relative difference
        is calculated by
            (x - y)/( (x + y)/2 ).
        x and y are equal when the relative difference is 0. When both x and y
        are 0, the mean ( (x + y)/2 ) is set to 1 to ensure the result is 0.
        
        Args:
            mat1 (ndarray): 
                Original Hi-C matrix.
            mat2 (ndarray):
                Modified Hi-C matrix.

        Return:
            diff (n x n array):
                The relative difference between mat1 and mat2
        """
    # calculate the average IF
    mean = (mat1 + mat2)/2

    # make sure that pairs where both IF are zero returns 0
    mean[mean == 0] = 1

    # calculate the relative difference
    diff = (mat1 - mat2)/mean
    
    return diff

def ratio(mat1, mat2):
    """Calculates the ratio of two Hi-C matrices.
        
        The calculations are done in the same way as in HiGlass' divide by.
        All division by zero is set to NaN.
        
        Args:
            mat1 (ndarray): 
                Original Hi-C matrix.
            mat2 (ndarray):
                Modified Hi-C matrix.

        Return:
            diff (n x n array):
                The ratio of mat1 and mat2.
        """
    # replaces all zeros in denominator to NaN to ensure no division by zero
    mat2[mat2 == 0] = np.nan

    diff = mat1/mat2

    # add pseudocount - the smallest value in diff
    #diff += np.nanmin(diff)
    """
    if replace_zero_zero:
        # find all 0/0
        sum_data = mat1 + mat2

        # set all 0/0 to 1
        diff[sum_data == 0] = 1
    """
    return diff