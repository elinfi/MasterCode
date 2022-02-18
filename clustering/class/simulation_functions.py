import re
import cooler

import numpy as np
import pyranges as pr

def concat_region(chromosome, start, end):
    """Concatenate region string.
    
    Args:
        chromosome (string):
            Chromosome name in the style chrX.
        start (int):
            Start position.
        end (int):
            End position.
    Returns:
        region (string):
            Genomic range string in the style {chrom}:{start}-{end}.
    """
    # concatenate region string
    region = chromosome + ':' + str(int(start)) + '-' + str(int(end))
    
    return region

def split_region(region):
    chrom, start, end = re.split("\W+", region)
    start, end = int(start), int(end)
    
    return chrom, start, end
    
    
def get_region(max_range):
    """Creates random region with range same as max_range.
    
    Args:
        max_range (int):
            Size of genomic range.
    Returns:
        region (string):
            Genomic range string in the style {chrom}:{start}-{end}.
    """
    # get random chromosome and corresponding chromosome size
    chrom, chr_size = random_chromosome()
    
    # get random start position
    start = np.random.randint(chr_size - max_range)
    end = start + max_range
    
    # concatenate genom range string
    region = concat_region(chrom, start, end)
    
    return region

def tad_df():
    """Extract tad dataframe from bed file.
    
    Args:
    
    Returns:
        df (pandas datafram):
            Dataframe containing TAD posititions.
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
    
    return df

def find_tads(region):
    """Find TADs within given region and store in dataframe.
    
    Args:
        region (string):
            Genomic range string in the style {chrom}:{start}-{end}.
            
    Returns:
        df (pandas df):
            Pandas dataframe containing tads within given region.
    """
    # split genomic range string
    chrom, start, end = split_region(region)
    
    # get TAD dataframe
    df = tad_df()
    
    # extract TADs for given chromosome
    df = df.loc[df['Chromosome'] == chrom]
    
    # extract TADs within region [start, end]
    df = df.loc[df['Start'] >= start]
    df = df.loc[df['End'] <= end]
    
    return df

def get_tad_region(region):
    """Choose random tad from df.
    
    Args:
        df (pandas dataframe):
            Pandas dataframe containing tads within a region.
    """
    df = find_tads(region)
    
    # choose random tad within given genomic range
    tad = df.sample(n=1)
    
    chrom = tad['Chromosome'].values[0]
    start = tad['Start'].values[0]
    end = tad['End'].values[0]
    
    region = concat_region(chrom, start, end)
    return region

def cooler_obj(filename, resolution):
    clr = cooler.Cooler(filename
                        + '::resolutions/'
                        + str(resolution))    
    return clr

def get_tad_idx(mat_region, tad_region, resolution):
    filename = '/home/elinfi/coolers/HiC_wt_merged.mcool'
    clr = cooler_obj(filename, resolution)
    
    mat_extent = clr.extent(mat_region)
    tad_extent = clr.extent(tad_region)
    
    tad_i = tad_extent[0] - mat_extent[0]
    tad_j = tad_extent[1] - mat_extent[1]
    
    return tad_i, tad_j

def matrix(region, replicate, resolution):
    filename = '/home/elinfi/coolers/HiC_' + replicate + '.mcool'
    clr = cooler_obj(filename, resolution)
    matrix = clr.matrix(balance=True).fetch(region)
    
    return matrix
    






































































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
    elif method == 'subtract':
        matrix = wt_002.subtract(wt_001)
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