"""Simulation of Hi-C data

Creates empty matrix with background noise.
"""
import cooler

import numpy as np
import pyranges as pr

from data_preparation import DataPreparation


def pick_tad(maxSize):
    """Get region for random TAD below a given size.
    
    Args:
        maxSize (int):
            The maximum size of TAD region
    Returns:
        region (string):
            Genomic range string in the style {chrom}:{start}-{end} containing
            TAD region.
    """
    # bed-file containing TAD positions for wild type merged
    filename = '/home/elinfi/storage/master_project/processed_data/tads/' \
               + 'HiC_wt_merged_normalized_and_corrected_ice_domains.bed'

    # read in the bed file with column names
    df = pr.read_bed(filename, as_df=True)

    # remove unecessary columns
    df = df.drop(columns=['Score', 'Strand', 'ThickStart', 'ThickEnd', 
                          'ItemRGB'])

    # add new column containing TAD sizes
    df['Size'] = df['End'] - df['Start']

    # Sort table by TAD size
    df = df.sort_values(['Size'])

    # get random TAD below a given size
    tad = df.loc[df['Size'] <= maxSize].sample(n=1)

    # create region
    region = tad['Chromosome'].values[0] \
             + ':' \
             + str(tad['Start'].values[0]) \
             + '-' \
             + str(tad['End'].values[0])

    return region

def get_tad(region, resolution):
    """Extract TAD from merged wild type at given resolution and region.
    
    Args:
        region (string):
            Genomic range string in the style {chrom}:{start}-{end} containing
            TAD region. Unit prefixes k, M, G are supported for start and end.
            
    Returns:
        matrix (ndarray):
            Matrix containing Hi-C data for given TAD.
    """
    filename = '/home/elinfi/coolers/HiC_wt_merged.mcool'
    clr = cooler.Cooler(filename
                        + '::resolutions/'
                        + str(resolution))
    matrix = clr.matrix(balance=True).fetch(region)
    
    return matrix

def df_chromsizes():
    """Extract chromosome sizes as pandas series.
    
    Args:
    
    Returns:
        df (pandas series):
            Dataframe containing all chromosome names and corresponding sizes.
    """
    filename = '/home/elinfi/coolers/HiC_wt_001_1000.cool'
    
    # get dataframe containing chromosome sizes
    df = cooler.Cooler(path_wt_001).chromsizes
    
    return df

def random_chromosome():
    """Extract random chromosome and corresponding chromosome size.
    
    Returns:
        chromosome (string):
            Chromosome name.
        size (int):
            Corresponding chromosome size.
    """
    # data frame with chromosome sizes
    df = df_chromsizes()
    
    # random index
    idx = np.random.randint(0, len(df))
    
    # get chromosome name and corresponding size
    chromosome = df.index[idx]
    size = df[idx]
    
    return chromosome, size

def noise_matrix(resolution, region):
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

    matrix = wt_002.higlass_ratio(wt_001)

    return matrix
    

def random_noise(resolution, max_range):
    """Create matrix with noise.
    
    Args:
        resolution (int):
            Resolution for cooler object.
        max_range (int):
            Size of genomic range.
    Returns:
        noise (ndarray):
            Matrix containing noise.
    """
    # get random chromosome and corresponding chromosome size
    chromosome, chr_size = random_chromosome()
    
    # get random start position
    start = np.random.randint(chr_size - max_range)
    end = start + max_range
    
    # concatenate region string
    region = chromosome + ':' + str(start) + '-' + str(end)
    
    # extract noise matrix for given region and resolution
    noise = noise_matrix(resolution, region)
    
    return noise

    

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
        self.resolution = resolution
        self.max_range = max_range
        
    
    
    