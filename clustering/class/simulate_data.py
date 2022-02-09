"""Simulation of Hi-C data

Creates empty matrix with background noise.
"""

import numpy as np
import pyranges as pr

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
    
    
    def extract_tads(self):
        # bed-file containing TAD positions for wild type merged
        filename = '/home/elinfi/storage/master_project/processed_data/tads/' \
                   + 'HiC_wt_merged_normalized_and_corrected_ice_domains.bed'
    
        
        df = pr.read_bed(filename, as_df=True)