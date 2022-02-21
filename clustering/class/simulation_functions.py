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

def get_region_with_ntads(max_range, ntads):
    """Check that random region of size max_range contains at least ntads TADs.
    
    Args:
        max_range (int):
            Size of genomic range.
        ntads (int):
            Minimum number of TADs to be within genomic range.
            
    Returns:
        region (string):
            Genomic range string in the style {chrom}:{start}-{end}.
    """
        
    test = True
    
    while test:
        region = get_region(max_range)
        df = find_tads(region)
        test = df.shape[0] < ntads
    
    return region, df
        
    
    
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

def df2tad(df, n=1, random_state=None):
    # choose random tad within given genomic range
    tad = df.sample(n=n, random_state=random_state)
    tad = tad.sort_values(['Start'])
       
    tad_regions = []
    for i in range(tad.shape[0]):
        chrom = tad['Chromosome'].values[i]
        start = tad['Start'].values[i]
        end = tad['End'].values[i]
        
        tad_region = concat_region(chrom, start, end)
        tad_regions.append(tad_region)
    return tad_regions

def region2tad(region):
    """Choose random tad from df.
    
    Args:
        df (pandas dataframe):
            Pandas dataframe containing tads within a region.
    """
    df = find_tads(region)
    
    tad_region = df2tad(df, n=1)
    
    return tad_region

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
    print(mat_extent)
    print(tad_extent)
    
    tad_i = tad_extent[0] - mat_extent[0]
    tad_j = tad_extent[1] - mat_extent[0]
    
    return tad_i, tad_j

def matrix(region, replicate, resolution):
    filename = '/home/elinfi/coolers/HiC_' + replicate + '.mcool'
    clr = cooler_obj(filename, resolution)
    matrix = clr.matrix(balance=True).fetch(region)
    
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