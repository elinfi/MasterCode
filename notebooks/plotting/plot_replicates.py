import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/class/')

import cooltools.lib.plotting

import seaborn as sns
import pretty_plotting as pplot
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, LogFormatterSciNotation

def plot_replicates(org, mod, region, 
                    s=18, m=24, l=26,
                    figsize=(16, 7.5),
                    background=False, 
                    save=False):
    if background:
        pplot.background_color('gray-lighter')
      
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # create subplots
    f, axs = plt.subplots(figsize=figsize,
                          nrows=1,
                          ncols=2,
                          sharex=True, sharey=True)
    
    # adjust subplots
    plt.subplots_adjust(left=0.1,
                        bottom=0.15, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.07)
    
    # norm and extent to used in plots
    norm = LogNorm()
    extent=pplot.region2extent(region)
    
    # wt 001
    ax = axs[0]
    im = ax.matshow(org,
                    cmap='fall',
                    norm=norm,
                    extent=extent)
    ax.set_title('WT 001', y=1.01)
    ax.tick_params(axis='x', labelsize=s)
    ax.tick_params(axis='y', labelsize=s)
    pplot.format_ticks(ax)
    
    # wt 002
    ax = axs[1]
    im = ax.matshow(mod,
                    cmap='fall',
                    norm=norm,
                    extent=extent)
    
    ax.set_title('WT 002', y=1.01)
    ax.tick_params(axis='x', labelsize=s)
    ax.tick_params(axis='y', labelsize=s)
    pplot.format_ticks(ax)
    
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction frequency', 
                 ax=list(axs), ticks=LogLocator(base=2), 
                 format=LogFormatterSciNotation(base=2))
    
    if save:
        plt.savefig(f'/home/elinfi/MasterCode/img/Replicates/{region}_wt_001_002.png', idp=500)
