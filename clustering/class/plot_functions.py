import os

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import TwoSlopeNorm

import pretty_plotting as pplot

def sub_colorbars(sub_data, region, filename, 
                  PATH='../../Figures/Differences/Colorbar', 
                  fontsizes=(22, 26, 30), figsize=(27,10),
                  sharex=True, sharey=False, hspace=0.4,
                  constrained_layout=True):
    """Plots comparision of subdata with different colorbars.
    
    The first colorbar is centered using the maximum value of the data as
    constraints, whereas the other colorbar is centered at 0 with vmin and vmax
    as the min and max value respectively.
    
    Args:
        sub_data (ndarray):
            Matrix containing subtraction difference Hi-C data.
        region (str):
            Genomic region string used in sub_data.
        filename (str):
            Name of output file.
        PATH (str):
            Path to storage of output file.
        fontsize (tuple):
            Fontsizes given on the form (small, medium, large).
        figsize (tuple):
            Size of figure.
        constrained_layout (boolean):
            Whether or not to use constrained layout in plotting.
    Returns:
        Saves figures in file given in outputname.
    """
    # set font sizes
    s, m, l = fontsizes
    pplot.font_size(s, m, l)
    
    # set axis extent based on genomic region
    extent = pplot.region2extent(region)
        
    # create figure
    f, axs = plt.subplots(figsize=figsize,
                          nrows=1,
                          ncols=2,
                          sharex=sharex, sharey=sharey,
                          squeeze=True,
                          constrained_layout=constrained_layout)
    #plt.subplots_adjust(hspace=hspace)
    plt.subplots_adjust(left=0.1,
                        bottom=0.15, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4)
                        #hspace=0.5)
    
    # first subplot
    ax = axs[0]
    vmax = np.nanmax(abs(sub_data))
    im = ax.matshow(sub_data,
                    vmax=vmax,
                    vmin=-vmax,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', ax=ax)
    ax.set(title='Subtraction')
    pplot.format_ticks(ax)
    
    # second subplot    
    ax = axs[1]
    im = ax.matshow(sub_data,
                    norm=TwoSlopeNorm(vcenter=0),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', ax=ax)
    ax.set(title='Subtraction')
    pplot.format_ticks(ax)
    
    # save figure
    plt.savefig(os.path.join(PATH, f'{region}_' + filename))
    
def sub_vmax(sub_data, region, filename, 
             PATH='../../Figures/Differences/Colorbar', 
             fontsizes=(16, 20, 24), figsize=(10, 8),
             sharex=True, sharey=False,
             left=0.1, bottom=0.15, right=0.9, top=0.9,
             constrained_layout=False):
    """Plot subtraction with vmax/vmin constraints on colorbar.
    """
    
    # set font sizes
    s, m, l = fontsizes
    pplot.font_size(s, m, l)
    
    # set axis extent based on genomic region
    extent = pplot.region2extent(region)
        
    # create figure
    f, axs = plt.subplots(figsize=figsize,
                          nrows=1,
                          ncols=1,
                          sharex=sharex, sharey=sharey,
                          squeeze=True,
                          constrained_layout=constrained_layout)
    #plt.subplots_adjust(hspace=hspace)
    plt.subplots_adjust(left=left,
                        bottom=bottom,
                        right=right,
                        top=top)
                        #hspace=0.5)
    
    # creates plot
    ax = axs
    vmax = np.nanmax(abs(sub_data))
    im = ax.matshow(sub_data,
                    vmax=vmax,
                    vmin=-vmax,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', ax=ax)
    ax.set(title='Subtraction')
    pplot.format_ticks(ax)
    
    # save figure
    plt.savefig(os.path.join(PATH, f'{region}_' + filename))
    
def sub_norm(sub_data, region, filename, 
             PATH='../../Figures/Differences/Colorbar', 
             fontsizes=(16, 20, 24), figsize=(10, 8),
             left=0.1, bottom=0.15, right=0.9, top=0.9,
             sharex=True, sharey=False,
             constrained_layout=False):
    """Plot subtraction with vmax/vmin constraints on colorbar.
    """
    
    # set font sizes
    s, m, l = fontsizes
    pplot.font_size(s, m, l)
    
    # set axis extent based on genomic region
    extent = pplot.region2extent(region)
        
    # create figure
    f, axs = plt.subplots(figsize=figsize,
                          nrows=1,
                          ncols=1,
                          sharex=sharex, sharey=sharey,
                          squeeze=True,
                          constrained_layout=constrained_layout)
    #plt.subplots_adjust(hspace=hspace)
    plt.subplots_adjust(left=left,
                        bottom=bottom,
                        right=right,
                        top=top)
                        #hspace=0.5)
    
    # creates plot
    ax = axs
    im = ax.matshow(sub_data,
                    norm=TwoSlopeNorm(vcenter=0),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', ax=ax)
    ax.set(title='Subtraction')
    pplot.format_ticks(ax)
    
    # save figure
    plt.savefig(os.path.join(PATH, f'{region}_' + filename))