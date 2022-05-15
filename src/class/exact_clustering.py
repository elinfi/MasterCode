"""Simulation of Hi-C data

Creates empty matrix with background noise.
"""
import cooler

import numpy as np
import comparison as comp
import pyranges as pr
import simulation_functions as sim

from data_preparation import DataPreparation


class ExactClustering():
    def __init__(self, mat_region='chr10:6351511-10351511', resolution=16000):
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
        self.mat_region = mat_region
        self.resolution = resolution
        
        self.tad_df = sim.find_tads(self.mat_region)
        self.loop_df = sim.find_loop(self.mat_region)
        
        org = sim.matrix(self.mat_region, 'wt_001', self.resolution)
        mod = sim.matrix(self.mat_region, 'wt_002', self.resolution)
        
        self.org = np.ones(org.shape)
        self.mod = np.ones(mod.shape)
        
    def print_tads_loops(self):
        print(f'#TADs: {self.tad_df.shape[0]}')
        print(f'#Loops: {self.loop_df.shape[0]}')
    
    def change_tad(self, change, nums, **kwargs):
        """Change TAD in mod
        
        Args:
            change (function):
                Function for how to change the TAD.
                Example of function for doubling IF in TAD: lambda x: 2*x
        """
        tad_regions = sim.df2tad(self.tad_df, nums)
        
        for tad_region in tad_regions:
            i, j = sim.get_tad_idx(self.mat_region, 
                                   tad_region,
                                   self.resolution)

            self.mod[i:j, i:j] = change(self.mod[i:j, i:j], **kwargs)
            print(i, j)

            
    def change_tad_tad(self, change, nums, **kwargs):
        """Change TAD-TAD interactions in mod.
        
        Args:
            change (function):
                Function for how to change the TAD.
                Example of function for doubling IF in TAD: lambda x: 2*x
            nums (list):
                List of two TAD indices.
        """
        tad_regions = sim.df2tad(self.tad_df, nums)
        
        i1, j1 = sim.get_tad_idx(self.mat_region, 
                                 tad_regions[0],
                                 self.resolution)
        i2, j2 = sim.get_tad_idx(self.mat_region, 
                                 tad_regions[1],
                                 self.resolution)

        self.mod[i1:j1, i2:j2] = change(self.mod[i1:j1, i2:j2], **kwargs)
        self.mod[i2:j2, i1:j1] = change(self.mod[i2:j2, i1:j1], **kwargs)
        
    def change_loop(self, change, nums, **kwargs):
        loops = sim.df2loop(self.loop_df, nums)
        for i in range(len(loops)):
            i1, j1 = sim.get_tad_idx(self.mat_region,
                                     loops[i][0],
                                     self.resolution)
            i2, j2 = sim.get_tad_idx(self.mat_region,
                                     loops[i][1],
                                     self.resolution)
            
            self.mod[i1:j1, i2:j2] = change(self.mod[i1:j1, i2:j2], **kwargs)
            self.mod[i2:j2, i1:j1] = change(self.mod[i2:j2, i1:j1], **kwargs)
            print(i1)
            
    def change_stripe(self, i, j, l, change, **kwargs):
        """
        Regulate stripe interaction.
        
        Args:
            i (int):
                left and top
        """
        depth = 2
        self.mod[i:j, l-depth:l] = change(self.mod[i:j, l-depth:l], **kwargs)
        self.mod[l-depth:l, i:j] = change(self.mod[l-depth:l, i:j], **kwargs)
        
    def compare(self, method, **kwargs):
        if method == 'ratio':
            diff = comp.ratio(self.org.copy(), self.mod.copy())
        elif method == 'reldiff':
            diff = comp.relative_difference(self.org.copy(), self.mod.copy())
        elif method == 'reldiffmax':
            diff = comp.relative_difference_max(self.org.copy(), self.mod.copy())
        elif method == 'reldiff_sqrt':
            diff = comp.relative_difference_sqrt(self.org.copy(), self.mod.copy())
        elif method == 'sub':
            diff = comp.subtraction(self.mod.copy(), self.org.copy())
        elif method == 'pseudo':
            diff = comp.pseudo_diff(self.mod.copy(), self.org.copy(), **kwargs)
        else:
            raise NameError("The method is not valid. Use 'ratio' or 'reldiff'")
        
        return diff
    
    def write2file(self, filename, matrix):
        np.save(filename, matrix)
        
    
        
        
    
    
    