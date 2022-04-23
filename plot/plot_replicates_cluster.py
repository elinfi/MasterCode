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

PATH_DIFF = '/home/elinfi/MasterCode/data/replicates/'
PATH_CLUSTER = '/home/elifni/MasterCode/data/clusters/'
REGION = 'chr10:6351511-10351511'

################################################################################
################################################################################
################################################################################

# use seaborn style
sns.set_theme('paper')
sns.set_style('ticks')

################################################################################
# FOLD CHANGE ##################################################################
################################################################################

diff = np.load(os.path.join(PATH, f'{REGION}_ratio.npy'))
cluster_lfc = np.load(os.path.join(PATH_CLUSTER, f'{REGION}_ratio_log_cluster_wdiag_0_medoids_5.npy'))
cluster_higlass = np.load(os.path.join(PATH_CLUSTER, f'{REGION}_ratio_higlass_cluster_wdiag_0_medoids_5.npy'))
    
# set font sizes
pplot.font_size(16, 18, 20)

# set extent and ticks label
extent = pplot.region2extent(region)    
locator2 = LogLocator(base=2)
formatter2 = LogFormatterSciNotation(base=2)

# create subplot
f, axs = plt.subplots(figsize=(16, 16),
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
ax = axs[0, 0]
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
ax = axs[0, 1]

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


ax = axs[1, 0]

if save:
    plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', path, 
                             f'{region}_fold_change_two.png'))