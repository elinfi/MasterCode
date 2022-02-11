"""Simulation of Hi-C data

Creates empty matrix with background noise.
"""
import cooler

import numpy as np
import pyranges as pr

from data_preparation import DataPreparation


def pick_tad(max_range):
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
    tad = df.loc[df['Size'] <= max_range].sample()

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
    df = cooler.Cooler(filename).chromsizes
    
    return df


def noise_matrix(resolution, region, method='ratio'):
    """
    Creates a matrix with background noise.

    Args:
        resolution (int):
            Resolution for cooler object
        region  (string): 
            Genomic range string of the style {chrom}:{start}-{end}, unit 
            prefixes k, M, G are supported
        method (string):
            ['ratio', 'reldiff']
            Method for calculating noise.

    Returns:
        matrix (ndarray):
            Matrix with the remaining noise difference from the two 
            replicates wild type 001 and wild type 002.
    """
    path_wt_001 = '/home/elinfi/coolers/HiC_wt_001.mcool'
    path_wt_002 = '/home/elinfi/coolers/HiC_wt_002.mcool'

    wt_001 = DataPreparation(path_wt_001, resolution, region)
    wt_002 = DataPreparation(path_wt_002, resolution, region)
    
    if method == 'ratio':
        matrix = wt_002.higlass_ratio(wt_001)
    elif method == 'reldiff':
        matrix = wt_002.relative_difference(wt_001)
    else:
        raise NameError("The method is not valid. Use 'ratio' or 'reldiff'")

    return matrix
    
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

def random_noise(resolution, max_range, method='ratio'):
    """Create matrix with noise.
    
    Args:
        resolution (int):
            Resolution for cooler object.
        max_range (int):
            Size of genomic range.
            method (string):
            ['ratio', 'reldiff']
            Method for calculating noise.
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
    noise = noise_matrix(resolution, region, method)
    
    return noise, region

def relative_difference(mat1, mat2):
    # calculate the average IF
    mean = (mat1 + mat2)/2

    # make sure that pairs where both IF are zero returns 0
    mean[mean == 0] = 1

    # calculate the relative difference
    diff = (mat1 - mat2)/mean
    return diff

def ratio(mat1, mat2):
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



class SimulateData():
    def __init__(self, resolution, max_range, method):
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
        self.max_range = int(max_range)
        self.method = method
        self.noise, self.region = random_noise(self.resolution, self.max_range, self.method)
        self.mat1 = np.zeros_like(self.noise)
        self.mat2 = np.zeros_like(self.noise)
        
    def add_tad(self):
        region = pick_tad(self.max_range)
        tad = get_tad(region, self.resolution)

        tad_size = tad.shape[0]
        max_idx = self.noise.shape[0] - tad_size

        # get random position for tad along the diagonal
        start = np.random.randint(low=0, high=max_idx) # [low, high)
        end = start + tad_size

        # add tad to simulation
        self.mat1[start:end, start:end] = tad
        self.mat2[start:end, start:end] = 4*tad
        
    def compare(self):
        if self.method == 'ratio':
            diff = ratio(self.mat1, self.mat2)
        elif self.method == 'reldiff':
            diff = relative_difference(self.mat1, self.mat2)
        else:
            raise NameError("The method is not valid. Use 'ratio' or 'reldiff'")
            
        final =self.noise + diff
        
        return final
        
    
        
        
    
    
    