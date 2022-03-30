import re 

import simulation_functions as sim
import matplotlib.pyplot as plt

from matplotlib.colors import to_rgba
from matplotlib.ticker import EngFormatter

def format_ticks(ax, x=True, y=True, rotate=True):
    bp_formatter = EngFormatter('b')
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)
    return None

def get_range(region):
    prefixes = {'k':1e3, 'M':1e6, 'G':1e9}
    
    try:
        num_region = int(region)
    except ValueError:
        value = re.match("\d+", region).group()
        prefix = re.search(r"[a-zA-Z]", region).group()
        num_region = int(float(value)*prefixes[prefix])
        
    return num_region
        

def split_region(region):
    chrom, start, end = re.split("\W+", region)
    
    start = get_range(start)
    end = get_range(end)
    
    return chrom, start, end

def region2extent(region):
    chrom, start, end = split_region(region)
    
    extent = (start, end, end, start)
    
    return extent

def background_color(color='gray-light'):
    colors = {'gray-light': '#cccccc', 
              'gray-lighter': '#e5e5e5', 
              'gray-lightest': '#f2f2f2'}
    
    plt.rcParams['axes.facecolor'] = colors[color]

def font_size(SMALL_SIZE=12, MEDIUM_SIZE=14, BIGGER_SIZE=16):
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
