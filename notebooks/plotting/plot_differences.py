import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/class/')

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pretty_plotting as pplot

from matplotlib.ticker import LogLocator, LogFormatterSciNotation
from mid_point_log_norm import MidPointLogNorm
from matplotlib.colors import TwoSlopeNorm, LogNorm

def plot_replicates(org, mod, region, save=False):
    f, axs = plt.subplots(figsize=(23, 10),
                          nrows=1,
                          ncols=2,
                          sharex=True, sharey=False)

    norm = mpl.colors.LogNorm()
    extent=pplot.region2extent(region)
    pplot.background_color()
    s = 18
    pplot.font_size(s, 24, 26)
    plt.subplots_adjust(left=0.1,
                        bottom=0.15, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4)

    ax = axs[0]
    im = ax.matshow(org,
                    cmap='fall',
                    norm=norm,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction frequency', ax=ax,
                 ticks=LogLocator(base=2), format=LogFormatterSciNotation(base=2))
    ax.set_title('WT 001', y=1.01)
    ax.tick_params(axis='x', labelsize=s)
    ax.tick_params(axis='y', labelsize=s)
    pplot.format_ticks(ax)
    pplot.background_color()

    ax = axs[1]
    im = ax.matshow(mod,
                    cmap='fall',
                    norm=norm,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction frequency', ax=ax,
                 ticks=LogLocator(base=2), format=LogFormatterSciNotation(base=2))
    ax.set_title('WT 002', y=1.01)
    ax.tick_params(axis='x', labelsize=s)
    ax.tick_params(axis='y', labelsize=s)
    pplot.format_ticks(ax)
    pplot.background_color()
    
    if save:
        plt.savefig('../img/Replicates/wt_001_002.png')

def plot_comparison(org, mod, region, save=False):
    f, axs = plt.subplots(figsize=(20, 20),
                          nrows=1,
                          ncols=2,
                          sharex=True, sharey=False)

    norm = mpl.colors.LogNorm()
    extent=pplot.region2extent(region)
    pplot.background_color()

    ax = axs[0]
    im = ax.matshow(org,
                    cmap='fall',
                    norm=norm,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Log2 ratio difference', ax=ax,
                 ticks=LogLocator(base=2), format=LogFormatterSciNotation(base=2))
    ax.set_title('Original: ' + region, y=1.01)
    pplot.format_ticks(ax)

    ax = axs[1]
    im = ax.matshow(mod,
                    cmap='fall',
                    norm=norm,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Log2 ratio difference', ax=ax,
                 ticks=LogLocator(base=2), format=LogFormatterSciNotation(base=2))
    ax.set_title('Modified: ' + region, y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig('../img/Simulation/org_mod.png')
    
def plot_diff(diff, method, region, save=False, **kwargs):
    if method == 'reldiff' or method == 'reldiffmax':
        f, axs = plt.subplots(figsize=(10, 8),
                              nrows=1,
                              ncols=1,
                              sharex=True, sharey=False,
                              constrained_layout=True)
        ax = axs
        vmax = np.nanmax(abs(diff))
        extent = pplot.region2extent(region)
        norm = TwoSlopeNorm(vcenter=0)
        
        pplot.font_size(14, 18, 20)

        im = ax.matshow(diff,
                        #vmax=vmax,
                        #vmin=-vmax,
                        norm = norm,
                        cmap='bwr',
                        extent=extent)
        plt.colorbar(im, fraction=0.046, pad=0.04, label='Interactions difference', ax=ax)
        ax.set_title("Relative difference", y=1.01)
        pplot.format_ticks(ax)
        pplot.background_color()
        if save:
            plt.savefig(f'../img/Differences/{region}_reldiff_norm.png')
            
    if method == 'reldiff_sqrt':
        f, axs = plt.subplots(figsize=(10, 8),
                              nrows=1,
                              ncols=1,
                              sharex=True, sharey=False,
                              constrained_layout=True)
        ax = axs
        vmax = np.nanmax(abs(diff))
        extent = pplot.region2extent(region)
        norm = TwoSlopeNorm(vcenter=0)
        
        pplot.font_size(14, 18, 20)

        im = ax.matshow(diff,
                        #vmax=vmax,
                        #vmin=-vmax,
                        norm = norm,
                        cmap='bwr',
                        extent=extent)
        plt.colorbar(im, fraction=0.046, pad=0.04, label='Interactions difference', ax=ax)
        ax.set_title("Relative difference sqrt", y=1.01)
        pplot.format_ticks(ax)
        pplot.background_color()
        if save:
            plt.savefig(f'../img/Differences/{region}_reldiff_sqrt_norm.png')
            
    elif method == 'ratio':
        pplot.font_size(16, 20, 24)
        pplot.background_color()
        
        extent = pplot.region2extent(region)
        
        pseudo = np.nanmin(diff[np.nonzero(diff)])
        print(pseudo)
        vmin = np.nanmin(diff)
        vmax = np.nanmax(diff)
        print(np.sign(vmin + pseudo), vmax)
        norm=MidPointLogNorm(vmin=vmin+pseudo, vmax=vmax+pseudo, midpoint=1 + pseudo)
        locator2 = LogLocator(base=2)
        formatter2 = LogFormatterSciNotation(base=2)

        f, axs = plt.subplots(figsize=(14, 8),
                              nrows=1,
                              ncols=1,
                              sharex=True, sharey=False)
        ax = axs
        im = ax.matshow(diff + pseudo,
                        norm=norm,
                        cmap='bwr',
                        extent=extent)
        plt.colorbar(im, pad=0.04, label='$\log_2$ transformation', ax=ax,
                     ticks=locator2, format=formatter2)
        ax.set_title('Division')
        pplot.format_ticks(ax)
        if save:
            plt.savefig(f'../img/Differences/{region}_div_logs.png')
    
    elif method == 'pseudo':
        pplot.font_size(16, 20, 24)
        pplot.background_color()
        
        extent = pplot.region2extent(region)
        
        norm=MidPointLogNorm(midpoint=1)
        locator2 = LogLocator(base=2)
        formatter2 = LogFormatterSciNotation(base=2)

        f, axs = plt.subplots(figsize=(14, 8),
                              nrows=1,
                              ncols=1,
                              sharex=True, sharey=False)
        ax = axs
        im = ax.matshow(diff,
                        norm=norm,
                        cmap='bwr',
                        extent=extent)
        plt.colorbar(im, pad=0.04, label='$Interaction difference', ax=ax)
                     #ticks=locator2, format=formatter2)
        ax.set_title(f'Pseudocount = 1')
        pplot.format_ticks(ax)
        if save:
            plt.savefig(f'../img/Differences/{region}_pseudocount.png')
    