import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/class/')

import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pretty_plotting as pplot

from matplotlib.ticker import LogLocator, LogFormatterSciNotation
from mid_point_log_norm import MidPointLogNorm
from matplotlib.colors import TwoSlopeNorm, LogNorm


def plot_comparison(org, mod, region, save=False):
    f, axs = plt.subplots(figsize=(20, 20),
                          nrows=1,
                          ncols=2,
                          sharex=True, sharey=False)

    norm = mpl.colors.LogNorm()
    extent=pplot.region2extent(region)

    ax = axs[0]
    im = ax.matshow(org,
                    cmap='fall',
                    norm=norm,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Log2 ratio difference', ax=ax,
                 ticks=LogLocator(base=2), format=LogFormatterSciNotation(base=2))
    ax.set_title('Original: ' + region, y=1.01)
    pplot.format_ticks(ax)
    pplot.background_color(ax)

    ax = axs[1]
    im = ax.matshow(mod,
                    cmap='fall',
                    norm=norm,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Log2 ratio difference', ax=ax,
                 ticks=LogLocator(base=2), format=LogFormatterSciNotation(base=2))
    ax.set_title('Modified: ' + region, y=1.01)
    pplot.format_ticks(ax)
    pplot.background_color(ax)
    
    if save:
        plt.savefig('../img/Simulation/org_mod.png')

def plot_subtraction_norm(diff, region,
                          s=16, m=18, l=20,
                          figsize=(8.5, 7.5),
                          path='Replicates',
                          save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    norm = TwoSlopeNorm(vcenter=0)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                              nrows=1,
                              ncols=1,
                              sharex=True, sharey=False)
    # adjust subplots
    plt.subplots_adjust(right=0.85)
    
    # plot difference
    ax= axs
    im = ax.matshow(diff,
                    norm = norm,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)
    ax.set_title("Subtraction", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_sub_norm.png'))
        
def plot_reldiff_norm(diff, region, 
                      s=16, m=18, l=20,
                      figsize=(8.5, 7.5),
                      path='Replicates',
                      save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    norm = TwoSlopeNorm(vcenter=0)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                              nrows=1,
                              ncols=1,
                              sharex=True, sharey=False)
    # adjust subplots
    plt.subplots_adjust(right=0.85)
    
    # plot difference
    ax= axs
    im = ax.matshow(diff,
                    norm = norm,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)
    ax.set_title("Relative difference", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_reldiff_norm.png'))
        
def plot_poisson(diff, region, 
                 s=16, m=18, l=20,
                 figsize=(8.5, 7.5),
                 path='Replicates',
                 save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    norm = TwoSlopeNorm(vcenter=0)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                              nrows=1,
                              ncols=1,
                              sharex=True, sharey=False)
    # adjust subplots
    plt.subplots_adjust(right=0.85)
    
    # plot difference
    ax= axs
    im = ax.matshow(diff,
                    norm = norm,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)
    ax.set_title("Poisson", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_poisson.png'))
    
def plot_fold_change(diff, region, 
                     s=16, m=18, l=20,
                     figsize=(8.5, 7.5),
                     path='Replicates',
                     save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    norm = MidPointLogNorm(midpoint=1)
    
    locator2 = LogLocator(base=2)
    formatter2 = LogFormatterSciNotation(base=2)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                              nrows=1,
                              ncols=1,
                              sharex=True, sharey=False)
    # adjust subplots
    plt.subplots_adjust(right=0.85)
    
    # plot difference
    ax= axs
    im = ax.matshow(diff,
                    norm = norm,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax, ticks=locator2, format=formatter2)
    ax.set_title("Fold change", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_fold_change.png'))
    

def plot_fold_change_higlass(diff, region, 
                             s=16, m=18, l=20,
                             figsize=(8, 7.5),
                             path='Replicates',
                             save=False):
    vmin = np.nanmin(diff[diff > 0])
    vmax = np.nanmax(diff) + vmin
    
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    norm = MidPointLogNorm(midpoint=1+vmin)
    
    locator2 = LogLocator(base=2)
    formatter2 = LogFormatterSciNotation(base=2)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                              nrows=1,
                              ncols=1,
                              sharex=True, sharey=False)
    
    # plot difference
    ax= axs
    im = ax.matshow(diff + vmin,
                    norm = norm,
                    cmap='bwr',
                    vmin = vmin,
                    vmax = vmax,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax, ticks=locator2, format=formatter2)
    ax.set_title("HiGlass fold change", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_fold_change_higlass.png'))
        
def plot_fold_change_two(diff, region, 
                         s=16, m=18, l=20,
                         figsize=(16, 7.5),
                         path='Replicates',
                         save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and ticks label
    extent = pplot.region2extent(region)    
    locator2 = LogLocator(base=2)
    formatter2 = LogFormatterSciNotation(base=2)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                          nrows=1,
                          ncols=2,
                          sharex=True, sharey=True)
    
    # adjust subplots
    plt.subplots_adjust(#left=0.1,
                        #bottom=0.15, 
                        right=0.95, 
                        #top=0.9, 
                        wspace=0.2)
    
    # plot fold change
    ax = axs[0]
    im = ax.matshow(diff,
                    norm = MidPointLogNorm(midpoint=1),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax, ticks=locator2, format=formatter2)
    ax.set_title("Fold change", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # plot higlass
    ax = axs[1]
    
    vmin = np.nanmin(diff[diff > 0])
    vmax = np.nanmax(diff) + vmin
    
    im = ax.matshow(diff + vmin,
                    norm = MidPointLogNorm(midpoint=1+vmin),
                    cmap='bwr',
                    vmin = vmin,
                    vmax = vmax,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax, ticks=locator2, format=formatter2)
    ax.set_title("HiGlass fold change", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_fold_change_two.png'))
        
        
def plot_pseudocount(p0, p1, p2, p3,
                     region, 
                     s=16, m=18, l=20,
                     figsize=(18, 15),
                     path='Replicates',
                     save=False):

    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    norm = MidPointLogNorm(midpoint=1)
    
    # set extent and ticks label
    extent = pplot.region2extent(region)    
    locator2 = LogLocator(base=2)
    formatter2 = LogFormatterSciNotation(base=2)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                          nrows=2,
                          ncols=2,
                          sharex=True, sharey=True)
    
    # adjust subplots
    plt.subplots_adjust(left=0.07,
                        bottom=0.05, 
                        right=0.9, 
                        top=0.95, 
                        wspace=0.3,
                        hspace=0.15)
    
    # no pseudocount
    ax = axs[0, 0]
    im = ax.matshow(p0,
                    norm = MidPointLogNorm(midpoint=1),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)#, ticks=locator2, format=formatter2)
    ax.set_title("Fold change", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # pseudocount = 0.01
    ax = axs[0, 1]
    im = ax.matshow(p1,
                    norm = MidPointLogNorm(midpoint=1),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)#, ticks=locator2, format=formatter2)
    ax.set_title("Pseudocount = 0.01", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # pseudocount = 0.1
    ax = axs[1, 0]
    im = ax.matshow(p2,
                    norm = MidPointLogNorm(midpoint=1),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)#, ticks=locator2, format=formatter2)
    ax.set_title("Pseudocount = 0.1", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # pseudocount = 1
    ax = axs[1, 1]
    im = ax.matshow(p3,
                    norm = MidPointLogNorm(midpoint=1),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)#, ticks=locator2, format=formatter2)
    ax.set_title("Pseudocount = 1", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')

    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_pseudocounts_log10.png'))
        
        
def plot_pseudocount_one_colorbar(p0, p1, p2, p3,
                                  region, 
                                  s=16, m=18, l=20,
                                  figsize=(16, 15),
                                  path='Replicates',
                                  save=False):
    
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    vmax = np.max((np.nanmax(p0), np.nanmax(p1), np.nanmax(p2), np.nanmax(p3)))
    vmin = np.min((np.nanmin(p0[p0>0]), np.nanmin(p1[p1>0]), np.nanmin(p2[p2>0]), np.nanmin(p3[p3>0])))
    norm = MidPointLogNorm(vmin=vmin, midpoint=1, vmax=vmax)
    
    # set extent and ticks label
    extent = pplot.region2extent(region)    
    locator2 = LogLocator(base=2)
    formatter2 = LogFormatterSciNotation(base=2)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                          nrows=2,
                          ncols=2,
                          sharex=True, sharey=True)
    
    # adjust subplots
    plt.subplots_adjust(#left=0.1,
                        #bottom=0.15, 
                        right=0.95, 
                        #top=0.9, 
                        wspace=0)
    
    # no pseudocount
    ax = axs[0, 0]
    im = ax.matshow(p0,
                    norm = norm,
                    cmap='bwr',
                    extent=extent)
    ax.set_title("Fold change", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # pseudocount = 0.01
    ax = axs[0, 1]
    im = ax.matshow(p1,
                    norm = norm,
                    cmap='bwr',
                    extent=extent)
    ax.set_title("Pseudocount = 0.01", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # pseudocount = 0.1
    ax = axs[1, 0]
    im = ax.matshow(p2,
                    norm = norm,
                    cmap='bwr',
                    extent=extent)
    ax.set_title("Pseudocount = 0.1", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # pseudocount = 1
    ax = axs[1, 1]
    im = ax.matshow(p3,
                    norm = norm,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=axs.ravel().tolist(), ticks=locator2, format=formatter2)
    ax.set_title("Pseudocount = 1", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')

    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_pseudocounts_one_colorbar.png'))
        
    
def plot_diff(diff, method, region, 
              s=16, m=18, l=20,
              save=False, **kwargs):
            
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
    