
"""Simulation of Hi-C data

Creates empty matrix with background noise.
"""
import cooler

import numpy as np
import comparison as cf
import pyranges as pr
import simulation_functions as sim

from data_preparation import DataPreparation


class SimulateData():
    def __init__(self, resolution, max_range):
        """
        Simulate Hi-C data.
        
        Args:
            resolution (int):
                Resolution for cooler object
            max_range (int):
                Size of genomic range.
        Attrs:
            matrix (ndarray):
                Numpy zeroes array.
        """
        clr = sim.cooler_obj(resolution)
        self.mat_region = sim.get_region(max_range)
        self.tad_region = sim.get_tad_region(self.mat_region)
        
        self.mat1 = sim.matrix(clr, self.mat_region)
        self.mat2 = self.mat1.copy()
        
        self.tad_i, self.tad_j = sim.get_tad_idx(clr, 
                                                 self.mat_region, 
                                                 self.tad_region)
    
    def change_tad(self, change):
        """Change TAD in mat2
        
        Args:
            change (function):
                Function for how to change the TAD.
                Example of function for doubling IF in TAD: lambda x: 2*x
        """
        i, j = (self.tad_i, self.tad_j)
        
        self.mat2[i:j, i:j] = change(self.mat2[i:j, i:j])
        
    def compare(self, method):
        if method == 'ratio':
            diff = cf.ratio(self.mat1, self.mat2)
        elif method == 'reldiff':
            diff = cf.relative_difference(self.mat1, self.mat2)
        else:
            raise NameError("The method is not valid. Use 'ratio' or 'reldiff'")
        
        return diff
        
    
        
        
    
    
    