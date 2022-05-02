import os
import sys
import cooltools.lib.plotting
sys.path.insert(1, '/home/elinfi/MasterCode/src/class/')
sys.path.insert(2, '/home/elinfi/MasterCode/plot/func')

import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pretty_plotting as pplot

from copy import copy
from matplotlib.ticker import LogLocator, LogFormatterSciNotation, LogFormatterMathtext
from mid_point_log_norm import MidPointLogNorm
from matplotlib.colors import TwoSlopeNorm, LogNorm

PATH_DIFF = '/home/elinfi/MasterCode/data/simulations/comparison'
PATH_CLUSTER = '/home/elinfi/MasterCode/data/simulations/cluster/tad_2_3_4_stripe_2'
PATH_IMG = '/home/elinfi/MasterCode/img/simulations/pdf/clusters/tad_2_3_4_stripe_2'
REGION = 'chr10:6351511-10351511'
#EXTENSION = '_tad_1.4_2_1.7'
#EXTENSION = '_k_2'
EXTENSION = '_tad_2_3_4_stripe_2'
WDIAG = 0
MEDOIDS = 4
RSTATE = 19

################################################################################
# STYLE SETTINGS ###############################################################
################################################################################

# use seaborn style
sns.set_theme('paper')
sns.set_style('ticks')

# set extent and ticks label
extent = pplot.region2extent(REGION)    
locator2 = LogLocator(base=2)
#formatter2 = LogFormatterSciNotation(base=2)
formatter2 = LogFormatterMathtext(base=2)

# set font sizes
pplot.font_size(16, 18, 20)

# set colormap with ignored nan-values
bwr = copy(plt.cm.bwr)
bwr.set_bad('#cccccc', 1.0)

# create discrete colormap for clusters
clstr = copy(plt.get_cmap('cool', MEDOIDS))
clstr.set_bad('#cccccc', 1.0)
clstr_tick_locs = (np.arange(MEDOIDS) + 0.5)*(MEDOIDS-1)/MEDOIDS
clstr_tick_label = np.arange(1, MEDOIDS + 1)

################################################################################
# POISSON ######################################################################
################################################################################

diff = np.load(os.path.join(PATH_DIFF, REGION + EXTENSION + '_poisson.npy'))
clusters = np.load(os.path.join(PATH_CLUSTER, 
                                REGION + EXTENSION 
                                + f'_poisson_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.npy'))

# create subplot
f, axs = plt.subplots(figsize=(14.5, 6.5),
                      nrows=1,
                      ncols=2,
                      sharex=True, sharey=True)
# adjust subplots
plt.subplots_adjust(left=0.06,
                    bottom=0.12, 
                    right=0.96, 
                    top=0.945, 
                    wspace=0.2,
                    hspace=0.8)

# plot difference
ax = axs[0]
im = ax.matshow(diff,
                vmin = -np.nanmax(abs(diff)),
                vmax = np.nanmax(abs(diff)),
                cmap = bwr,
                extent = extent)
plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
             ax=ax)
ax.set_title("Poisson", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

# plot clusters
ax = axs[1]
im = ax.matshow(clusters,
                cmap=clstr,
                extent=extent)
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
cbar.set_ticks(clstr_tick_locs)
cbar.set_ticklabels(clstr_tick_label)

ax.set_title("Clustering", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

plt.savefig(os.path.join(PATH_IMG,
                         REGION + EXTENSION 
                         + f'_poisson_sym_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.pdf'))

################################################################################
# SUBTRACTION ##################################################################
################################################################################

diff = np.load(os.path.join(PATH_DIFF, REGION + EXTENSION + '_sub.npy'))
clusters = np.load(os.path.join(PATH_CLUSTER, 
                                REGION + EXTENSION 
                                + f'_sub_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.npy'))

# create subplot
f, axs = plt.subplots(figsize=(14.5, 6.5),
                      nrows=1,
                      ncols=2,
                      sharex=True, sharey=True)
# adjust subplots
plt.subplots_adjust(left=0.06,
                    bottom=0.12, 
                    right=0.96, 
                    top=0.945, 
                    wspace=0.2,
                    hspace=0.8)

# plot difference
ax = axs[0]
im = ax.matshow(diff,
                vmin = -np.nanmax(abs(diff)),
                vmax = np.nanmax(abs(diff)),
                cmap = bwr,
                extent = extent)
plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
             ax=ax)
ax.set_title("Subtraction", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

# plot clusters
ax = axs[1]
im = ax.matshow(clusters,
                cmap=clstr,
                extent=extent)
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
cbar.set_ticks(clstr_tick_locs)
cbar.set_ticklabels(clstr_tick_label)
ax.set_title("Clustering", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

plt.savefig(os.path.join(PATH_IMG,
                         REGION + EXTENSION 
                         + f'_sub_sym_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.pdf'))

################################################################################
# RELATIVE DIFFERENCE ##########################################################
################################################################################

diff = np.load(os.path.join(PATH_DIFF, REGION + EXTENSION + '_reldiff.npy'))
clusters = np.load(os.path.join(PATH_CLUSTER, 
                                REGION + EXTENSION 
                                + f'_reldiff_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.npy'))

# create subplot
f, axs = plt.subplots(figsize=(14.5, 6.5),
                      nrows=1,
                      ncols=2,
                      sharex=True, sharey=True)
# adjust subplots
plt.subplots_adjust(left=0.06,
                    bottom=0.12, 
                    right=0.96, 
                    top=0.945, 
                    wspace=0.2,
                    hspace=0.8)

# plot difference
ax = axs[0]
im = ax.matshow(diff,
                vmin = -np.nanmax(abs(diff)),
                vmax = np.nanmax(abs(diff)),
                cmap = bwr,
                extent = extent)
plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
             ax=ax)
ax.set_title("Relative difference", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

# plot clusters
ax = axs[1]
im = ax.matshow(clusters,
                cmap=clstr,
                extent=extent)
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
cbar.set_ticks(clstr_tick_locs)
cbar.set_ticklabels(clstr_tick_label)
ax.set_title("Clustering", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

plt.savefig(os.path.join(PATH_IMG,
                         REGION + EXTENSION 
                         + f'_reldiff_sym_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.pdf'))

################################################################################
# POISSON ######################################################################
################################################################################

c0 = np.load(os.path.join(PATH_CLUSTER, 
                          REGION + EXTENSION 
                          + f'_p{0}_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.npy'))
c1 = np.load(os.path.join(PATH_CLUSTER, 
                          REGION + EXTENSION 
                          + f'_p_{0.0001}_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.npy'))
c2 = np.load(os.path.join(PATH_CLUSTER, 
                          REGION + EXTENSION 
                          + f'_p_{0.001}_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.npy'))
c3 = np.load(os.path.join(PATH_CLUSTER, 
                          REGION + EXTENSION 
                          + f'_p_{1}_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.npy'))

# create subplot
f, axs = plt.subplots(figsize=(14.5, 13),
                      nrows=2,
                      ncols=2,
                      sharex=True, sharey=True)

# adjust subplots
plt.subplots_adjust(left=0.065,
                    bottom=0.065, 
                    right=0.96, 
                    top=0.97, 
                    wspace=0.15,
                    hspace=0.1)

# p = 0
ax = axs[0, 0]
im = ax.matshow(c0,
                cmap=clstr,
                extent=extent)
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
cbar.set_ticks(clstr_tick_locs)
cbar.set_ticklabels(clstr_tick_label)
ax.set_title("$p = 0$", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

# p = 0.0001
ax = axs[0, 1]
im = ax.matshow(c1,
                cmap=clstr,
                extent=extent)
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
cbar.set_ticks(clstr_tick_locs)
cbar.set_ticklabels(clstr_tick_label)
ax.set_title("$p = 0.0001$", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

# p = 0.001
ax = axs[1, 0]
im = ax.matshow(c2,
                cmap=clstr,
                extent=extent)
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
cbar.set_ticks(clstr_tick_locs)
cbar.set_ticklabels(clstr_tick_label)
ax.set_title("$p = 0.001$", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

# p = 1
ax = axs[1, 1]
im = ax.matshow(c3,
                cmap=clstr,
                extent=extent)
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, label='Clusters', ax=ax)
cbar.set_ticks(clstr_tick_locs)
cbar.set_ticklabels(clstr_tick_label)
ax.set_title("$p = 1$", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

plt.savefig(os.path.join(PATH_IMG,
                         REGION + EXTENSION 
                         + f'_pseudo_cluster_wdiag_{WDIAG}_medoids_{MEDOIDS}_rstate_{RSTATE}.pdf'))