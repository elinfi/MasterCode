import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/class/')

import cooltools.lib.plotting

import numpy as np
import seaborn as sns
import pretty_plotting as pplot
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, LogFormatterSciNotation
from mid_point_log_norm import MidPointLogNorm
from matplotlib.colors import TwoSlopeNorm, LogNorm


def sub_cluster_2(sub, cluster, medoids, region, wdiag,
                s=16, m=18, l=20,
                figsize=(18, 7.5),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
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
                        wspace=0.1)
    
    # plot difference
    ax = axs[0]
    im = ax.matshow(sub,
                    norm = TwoSlopeNorm(vcenter=0),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)
    ax.set_title("Subtraction", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # plot clusters
    ax = axs[1]
    # create discrete colormap
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    im = ax.matshow(cluster,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title("Clustering", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_sub_norm_wdiag_{wdiag}.png'))
        
def rel_cluster_2(diff, cluster, medoids, region, wdiag,
                s=16, m=18, l=20,
                figsize=(18, 7.5),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
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
                        wspace=0.1)
    
    # plot difference
    ax = axs[0]
    im = ax.matshow(diff,
                    norm = TwoSlopeNorm(vcenter=0),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)
    ax.set_title("Relative difference", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # plot clusters
    ax = axs[1]
    # create discrete colormap
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    im = ax.matshow(cluster,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title("Clustering", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_rel_norm_wdiag_{wdiag}.png'))
        
        
def poisson_cluster_2(diff, cluster, medoids, region, wdiag,
                s=16, m=18, l=20,
                figsize=(18, 7.5),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
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
                        wspace=0.1)
    
    # plot difference
    ax = axs[0]
    im = ax.matshow(diff,
                    norm = TwoSlopeNorm(vcenter=0),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)
    ax.set_title("Poisson", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # plot clusters
    ax = axs[1]
    # create discrete colormap
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    im = ax.matshow(cluster,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title("Clustering", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_poisson_norm_wdiag_{wdiag}.png'))
        
        
def ratio_cluster_2(diff, cluster, medoids, region, wdiag,
                s=16, m=18, l=20,
                figsize=(18, 7.5),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
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
                        wspace=0.1)
    
    # plot difference
    norm = MidPointLogNorm(midpoint=1)
    
    locator2 = LogLocator(base=2)
    formatter2 = LogFormatterSciNotation(base=2)
    
    ax = axs[0]
    im = ax.matshow(diff,
                    norm =norm,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax, ticks=locator2, format=formatter2)
    ax.set_title("Fold change", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # plot clusters
    ax = axs[1]
    # create discrete colormap
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    im = ax.matshow(cluster,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title("Clustering", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_ratio_log_wdiag_{wdiag}.png'))
        
        
def ratio_higlass_cluster_2(diff, cluster, medoids, region, wdiag,
                s=16, m=18, l=20,
                figsize=(18, 7.5),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
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
                        wspace=0.1)
    
    # plot difference
    vmin = np.nanmin(diff[diff > 0])
    vmax = np.nanmax(diff) + vmin
    norm = MidPointLogNorm(midpoint=1+vmin)
    
    locator2 = LogLocator(base=2)
    formatter2 = LogFormatterSciNotation(base=2)
    
    ax = axs[0]
    im = ax.matshow(diff + vmin,
                    norm =norm,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax, ticks=locator2, format=formatter2)
    ax.set_title("Fold change HiGlass", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # plot clusters
    ax = axs[1]
    # create discrete colormap
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    im = ax.matshow(cluster,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title("Clustering", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_ratio_higlass_wdiag_{wdiag}.png'))

################################################################################
################################################################################

def sub_cluster(sub, cluster1, cluster2, cluster3, medoids, region, wdiags,
                s=16, m=18, l=20,
                figsize=(18, 15),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                          nrows=2,
                          ncols=2,
                          sharex=True, sharey=True)
    # adjust subplots
    plt.subplots_adjust(left=0.07,
                        bottom=0.06, 
                        right=0.95, 
                        top=0.95, 
                        wspace=0.3,
                        hspace=0.15)
    
    # plot difference
    ax = axs[0, 0]
    im = ax.matshow(sub,
                    norm = TwoSlopeNorm(vcenter=0),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)
    ax.set_title("Subtraction", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    

    # create discrete colormap for cluster plotting
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax = axs[1, 0]
    
    # plot clusters
    ax = axs[0, 1]
    im = ax.matshow(cluster1,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[0]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 0]
    im = ax.matshow(cluster2,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[1]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 1]
    im = ax.matshow(cluster3,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[2]}", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_sub_norm_{wdiags}.png'))
        
def rel_cluster(diff, cluster1, cluster2, cluster3, medoids, region, wdiags,
                s=16, m=18, l=20,
                figsize=(18, 15),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                          nrows=2,
                          ncols=2,
                          sharex=True, sharey=True)
    # adjust subplots
    plt.subplots_adjust(left=0.07,
                        bottom=0.06, 
                        right=0.95, 
                        top=0.95, 
                        wspace=0.3,
                        hspace=0.15)
    
    # plot difference
    ax = axs[0, 0]
    im = ax.matshow(diff,
                    norm = TwoSlopeNorm(vcenter=0),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)
    ax.set_title("Relative difference", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # create discrete colormap for cluster plotting
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax = axs[1, 0]
    
    # plot clusters
    ax = axs[0, 1]
    im = ax.matshow(cluster1,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[0]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 0]
    im = ax.matshow(cluster2,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[1]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 1]
    im = ax.matshow(cluster3,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[2]}", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_rel_norm_{wdiags}.png'))
        
        
def poisson_cluster(diff, cluster1, cluster2, cluster3, medoids, region, wdiags,
                s=16, m=18, l=20,
                figsize=(18, 15),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                          nrows=2,
                          ncols=2,
                          sharex=True, sharey=True)
    # adjust subplots
    plt.subplots_adjust(left=0.07,
                        bottom=0.06, 
                        right=0.95, 
                        top=0.95, 
                        wspace=0.3,
                        hspace=0.15)
    
    # plot difference
    ax = axs[0, 0]
    im = ax.matshow(diff,
                    norm = TwoSlopeNorm(vcenter=0),
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)
    ax.set_title("Poisson", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # create discrete colormap for cluster plotting
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax = axs[1, 0]
    
    # plot clusters
    ax = axs[0, 1]
    im = ax.matshow(cluster1,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[0]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 0]
    im = ax.matshow(cluster2,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[1]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 1]
    im = ax.matshow(cluster3,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[2]}", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_poisson_norm_{wdiags}.png'))
        
        
def ratio_cluster(diff, cluster1, cluster2, cluster3, medoids, region, wdiags,
                s=16, m=18, l=20,
                figsize=(18, 15),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                          nrows=2,
                          ncols=2,
                          sharex=True, sharey=True)
    # adjust subplots
    plt.subplots_adjust(left=0.07,
                        bottom=0.06, 
                        right=0.95, 
                        top=0.95, 
                        wspace=0.3,
                        hspace=0.15)
    
    # plot difference
    norm = MidPointLogNorm(midpoint=1)
    
    locator2 = LogLocator(base=2)
    formatter2 = LogFormatterSciNotation(base=2)
    
    ax = axs[0, 0]
    im = ax.matshow(diff,
                    norm =norm,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax, ticks=locator2, format=formatter2)
    ax.set_title("Fold change", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # create discrete colormap for cluster plotting
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax = axs[1, 0]
    
    # plot clusters
    ax = axs[0, 1]
    im = ax.matshow(cluster1,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[0]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 0]
    im = ax.matshow(cluster2,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[1]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 1]
    im = ax.matshow(cluster3,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[2]}", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_ratio_log_{wdiags}.png'))
        
        
def ratio_higlass_cluster(diff, cluster1, cluster2, cluster3, medoids, region, wdiags,
                s=16, m=18, l=20,
                figsize=(18, 15),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                          nrows=2,
                          ncols=2,
                          sharex=True, sharey=True)
    # adjust subplots
    plt.subplots_adjust(left=0.07,
                        bottom=0.06, 
                        right=0.95, 
                        top=0.95, 
                        wspace=0.3,
                        hspace=0.15)
    
    # plot difference
    vmin = np.nanmin(diff[diff > 0])
    vmax = np.nanmax(diff) + vmin
    norm = MidPointLogNorm(midpoint=1+vmin)
    
    locator2 = LogLocator(base=2)
    formatter2 = LogFormatterSciNotation(base=2)
    
    ax = axs[0, 0]
    im = ax.matshow(diff + vmin,
                    norm =norm,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax, ticks=locator2, format=formatter2)
    ax.set_title("Fold change HiGlass", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # create discrete colormap for cluster plotting
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax = axs[1, 0]
    
    # plot clusters
    ax = axs[0, 1]
    im = ax.matshow(cluster1,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[0]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 0]
    im = ax.matshow(cluster2,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[1]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 1]
    im = ax.matshow(cluster3,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[2]}", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_ratio_higlass_{wdiags}.png'))
        
def pseudocount_cluster(diff, cluster1, cluster2, cluster3, medoids, region, wdiags, pseudo,
                s=16, m=18, l=20,
                figsize=(18, 15),
                path='Replicates/clusters',
                save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    # set extent and norm of plot
    extent = pplot.region2extent(region)
    
    # create subplot
    f, axs = plt.subplots(figsize=figsize,
                          nrows=2,
                          ncols=2,
                          sharex=True, sharey=True)
    # adjust subplots
    plt.subplots_adjust(left=0.07,
                        bottom=0.06, 
                        right=0.95, 
                        top=0.95, 
                        wspace=0.3,
                        hspace=0.15)
    
    # plot difference
    norm = MidPointLogNorm(midpoint=1)
    
    locator2 = LogLocator(base=2)
    formatter2 = LogFormatterSciNotation(base=2)
    
    ax = axs[0, 0]
    im = ax.matshow(diff,
                    norm =norm,
                    cmap='bwr',
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
                 ax=ax)#, ticks=locator2, format=formatter2)
    ax.set_title(f"Pseudocount = {pseudo}", y=1.01)
    
    pplot.format_ticks(ax)
    pplot.background_color(ax, 'gray-light')
    
    # create discrete colormap for cluster plotting
    cmap = plt.cm.cool
    bounds = np.linspace(0, medoids, medoids+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax = axs[1, 0]
    
    # plot clusters
    ax = axs[0, 1]
    im = ax.matshow(cluster1,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[0]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 0]
    im = ax.matshow(cluster2,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[1]}", y=1.01)
    pplot.format_ticks(ax)
    
    ax = axs[1, 1]
    im = ax.matshow(cluster3,
                    norm=norm,
                    cmap=cmap,
                    alpha=1,
                    extent=extent)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
    ax.set_title(f"w = {wdiags[2]}", y=1.01)
    pplot.format_ticks(ax)
    
    if save:
        plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                                 f'{region}_cluster_pseudo_{pseudo}_{wdiags}.png'))