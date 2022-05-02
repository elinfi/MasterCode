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

PATH_DIFF = '/home/elinfi/MasterCode/data/replicates/comparison'
PATH_REPS = '/home/elinfi/MasterCode/data/replicates'
PATH_IMG = 'replicates/pdf/comparison'
REGION = 'chr10:6351511-10351511'

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

fall = copy(plt.get_cmap('fall'))

################################################################################

wt1 = np.load(os.path.join(PATH_REPS, f'{REGION}_wt1.npy'))
wt2 = np.load(os.path.join(PATH_REPS, f'{REGION}_wt2.npy'))
wt1_zero = wt1 == 0
wt1_pos = wt1 > 0
wt2_zero = wt2 == 0
wt2_pos = wt2 > 0

wt1_zero = np.where(wt2_pos, wt1_zero, False)
wt2_zero = np.where(wt1_pos, wt2_zero, False)
print(np.sum(wt2_zero))

pseudo_raw = np.load(os.path.join(PATH_DIFF, f'{REGION}_p_0.0001.npy'))
pseudo = np.log2(pseudo_raw)

wt1_mean = np.nanmean(pseudo[wt2_zero])
wt2_mean = np.nanmean(pseudo[wt1_zero])
print(wt1_mean)
print(wt2_mean)

wt1_plot = np.where(~wt2_zero, np.nan, pseudo_raw)
wt2_plot = np.where(~wt1_zero, np.nan, pseudo_raw)
print(wt1_plot.shape)

# create subplots
f, axs = plt.subplots(figsize=(15, 7.5),
                      nrows=1,
                      ncols=2,
                      sharex=True, sharey=True)

# adjust subplots
plt.subplots_adjust(left=0.06,
                    bottom=0.06, 
                    right=0.94, 
                    top=0.995, 
                    wspace=0.04)

# wt 001
ax = axs[0]
vmin = min(np.nanmin(pseudo_raw), 1/np.nanmax(pseudo_raw))
vmax = max(np.nanmax(pseudo_raw), 1/np.nanmin(pseudo_raw))
im = ax.matshow(wt1_plot,
                cmap=bwr,
                norm = MidPointLogNorm(midpoint=1, vmin=vmin, vmax=vmax, log=np.log2),
                extent=extent)
ax.set_title('WT 001', y=1.01)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
pplot.format_ticks(ax)

# wt 002
ax = axs[1]
im = ax.matshow(wt2_plot,
                cmap=bwr,
                norm = MidPointLogNorm(midpoint=1, vmin=vmin, vmax=vmax, log=np.log2),
                extent=extent)

ax.set_title('WT 002', y=1.01)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
pplot.format_ticks(ax)

plt.colorbar(im, fraction=0.0235, pad=0.02, label='Interaction frequency', 
             ax=list(axs), ticks=LogLocator(base=2), 
             format=LogFormatterSciNotation(base=2))
    
plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', PATH_IMG, 
                         f'{REGION}_low_IF_region.pdf'))
