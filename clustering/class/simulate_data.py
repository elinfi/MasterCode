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
    def __init__(self, resolution, max_range, ntads):
        """
        Simulate Hi-C data.
        
        Args:
            resolution (int):
                Resolution for cooler object.
            max_range (int):
                Size of genomic range.
            ntads (int):
                Minimum number of TADs to be within genomic range.
        Attrs:
            matrix (ndarray):
                Numpy zeroes array.
        """
        self.resolution = resolution
        self.mat_region, self.df = sim.get_region_with_ntads(max_range, ntads)
        #self.mat_region = 'chr8:85013332-95013332'
        #self.tad_region = sim.get_tad_region(self.mat_region)
        #self.tad_region = 'chr8:86550000-87600000'
        
        self.mat1 = sim.matrix(self.mat_region, 'wt_001', self.resolution)
        self.mat2 = sim.matrix(self.mat_region, 'wt_002', self.resolution)
    
    def change_tad(self, change, n=1, random_state=None, **kwargs):
        """Change TAD in mat2
        
        Args:
            change (function):
                Function for how to change the TAD.
                Example of function for doubling IF in TAD: lambda x: 2*x
        """
        tad_regions = sim.df2tad(self.df, 1, random_state)
        
        for tad_region in tad_regions:
            i, j = sim.get_tad_idx(self.mat_region, 
                                   tad_region,
                                   self.resolution)

            self.mat2[i:j, i:j] = change(self.mat2[i:j, i:j], **kwargs)

            
    def change_tad_tad(self, change, random_state=None, **kwargs):
        tad_regions = sim.df2tad(self.df, 2, random_state)
        
        i1, j1 = sim.get_tad_idx(self.mat_region, 
                                 tad_regions[0],
                                 self.resolution)
        i2, j2 = sim.get_tad_idx(self.mat_region, 
                                 tad_regions[1],
                                 self.resolution)

        self.mat2[i1:j1, i2:j2] = change(self.mat2[i1:j1, i2:j2], **kwargs)
        self.mat2[i2:j2, i1:j1] = change(self.mat2[i2:j2, i1:j1], **kwargs)
        
    def compare(self, method):
        if method == 'ratio':
            diff = cf.ratio(self.mat1, self.mat2)
        elif method == 'reldiff':
            diff = cf.relative_difference(self.mat1, self.mat2)
        else:
            raise NameError("The method is not valid. Use 'ratio' or 'reldiff'")
        
        return diff
    
    def write2file(self, filename, matrix):
        np.save(filename, matrix)
        
    
        
        
    
    
    