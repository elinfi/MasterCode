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
fall.set_bad('w', 1.0)

################################################################################
# REPLICATES ###################################################################
################################################################################

wt1 = np.load(os.path.join(PATH_REPS, f'{REGION}_wt1.npy'))
wt2 = np.load(os.path.join(PATH_REPS, f'{REGION}_wt2.npy'))

print(np.nanmin(wt1[wt1 > 0]), np.nanmin(wt2[wt2 > 0]))
print(np.nanmax(wt1), np.nanmax(wt2))
print(np.nanmedian(wt1), np.nanmedian(wt2))
print(np.median(wt1[wt1 > 0]), np.median(wt2[wt2 > 0]))
boolean = wt2 == np.nanmin(wt2[wt2 > 0])
boolean = wt2 == 0
print(wt1[boolean].shape)
min_test = np.nanmin(wt1[boolean][wt1[boolean] > 0])
min2 = wt2[boolean][0]
min1 = wt1[boolean][0]
p = 0.0001
print(np.log2((min_test + p)/(0 + p)))


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
im = ax.matshow(wt1,
                cmap=fall,
                norm= LogNorm(),
                extent=extent)
ax.set_title('WT 001', y=1.01)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
pplot.format_ticks(ax)

# wt 002
ax = axs[1]
im = ax.matshow(wt2,
                cmap=fall,
                norm=LogNorm(),
                extent=extent)

ax.set_title('WT 002', y=1.01)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
pplot.format_ticks(ax)

plt.colorbar(im, fraction=0.0235, pad=0.02, label='Interaction frequency', 
             ax=list(axs), ticks=LogLocator(base=2), 
             format=LogFormatterSciNotation(base=2))
    
plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', PATH_IMG, 
                         f'{REGION}_wt1_wt2.pdf'))

################################################################################
# FOLD CHANGE ##################################################################
################################################################################

diff = np.load(os.path.join(PATH_DIFF, f'{REGION}_ratio.npy'))

# create subplot
fig, axs = plt.subplots(figsize=(16.8, 7.5),
                      nrows=1,
                      ncols=2,
                      sharex=True, sharey=True)

# adjust subplots
plt.subplots_adjust(left=0.04,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.95, 
                    wspace=0.1)

# plot fold change
ax = axs[0]
vmin = min(np.nanmin(diff[diff > 0]), 1/np.nanmax(diff))
vmax = max(np.nanmax(diff), 1/np.nanmin(diff[diff > 0]))
im = ax.matshow(diff,
                norm = MidPointLogNorm(midpoint=1, vmin=vmin, vmax=vmax, log=np.log2),
                cmap = bwr,
                extent=extent)
plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
             ax=ax, ticks=locator2, format=formatter2)
ax.set_title("Log fold change", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

# plot higlass
ax = axs[1]

vmin = np.nanmin(diff[diff > 0])
vmax = np.nanmax(diff) + vmin

im = ax.matshow(diff + vmin,
                norm = MidPointLogNorm(midpoint=1+vmin),
                cmap = bwr,
                vmin = vmin,
                vmax = vmax,
                extent = extent)
plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
             ax=ax, ticks=locator2, format=formatter2)
ax.set_title("HiGlass divide by", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', PATH_IMG, 
                         f'{REGION}_fold_change_two_sym.pdf'))

################################################################################
# RELATIVE DIFFERENCE ##########################################################
################################################################################

diff = np.load(os.path.join(PATH_DIFF, f'{REGION}_reldiff.npy'))

# create subplot
f, axs = plt.subplots(figsize=(8.8, 7.5),
                          nrows=1,
                          ncols=1,
                          sharex=True, sharey=False)

# adjust subplots
plt.subplots_adjust(left=0.01,
                    bottom=0.11, 
                    right=0.89, 
                    top=0.94)

# plot difference
ax= axs
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

plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', PATH_IMG, 
                         f'{REGION}_reldiff_sym.pdf'))

################################################################################
# SUBTRACTION ##################################################################
################################################################################

diff = np.load(os.path.join(PATH_DIFF, f'{REGION}_sub.npy'))

# create subplot
f, axs = plt.subplots(figsize=(9, 7.5),
                          nrows=1,
                          ncols=1,
                          sharex=True, sharey=False)

# adjust subplots
plt.subplots_adjust(left=0.1,
                    bottom=0.11, 
                    right=0.86, 
                    top=0.94)

# plot difference
ax= axs
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

plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', PATH_IMG,
                         f'{REGION}_sub_sym.pdf'))

################################################################################
# POISSON ######################################################################
################################################################################

diff = np.load(os.path.join(PATH_DIFF, f'{REGION}_poisson.npy'))

# create subplot
f, axs = plt.subplots(figsize=(9, 7.5),
                          nrows=1,
                          ncols=1,
                          sharex=True, sharey=False)

# adjust subplots
plt.subplots_adjust(left=0.1,
                    bottom=0.11, 
                    right=0.86, 
                    top=0.94)

# plot difference
ax = axs
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

plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', PATH_IMG,
                         f'{REGION}_poisson_sym.pdf'))

################################################################################
# PSEUDOCOUNT 0, 0.0001, 0.001, 1 ##############################################
################################################################################

p0 = np.load(os.path.join(PATH_DIFF, f'{REGION}_ratio.npy'))
p1 = np.load(os.path.join(PATH_DIFF, f'{REGION}_p_0.0001.npy'))
p2 = np.load(os.path.join(PATH_DIFF, f'{REGION}_p_0.001.npy'))
p3 = np.load(os.path.join(PATH_DIFF, f'{REGION}_p_1.npy'))

print(np.sum(np.isnan(p0)))
print(np.sum(np.isnan(p1)))
print(np.sum(np.isnan(p2)))
print(np.sum(np.isnan(p3)))

# create subplot
f, axs = plt.subplots(figsize=(14.2, 11.9),
                      nrows=2,
                      ncols=2,
                      sharex=True, sharey=True)

# adjust subplots
plt.subplots_adjust(left=0.05,
                    bottom=0.065, 
                    right=0.93, 
                    top=0.97, 
                    wspace=0.15,
                    hspace=0.1)

# no pseudocount
ax = axs[0, 0]
vmin = min(np.nanmin(p0[p0 > 0]), 1/np.nanmax(p0))
vmax = max(np.nanmax(p0), 1/np.nanmin(p0[p0 > 0]))
print(vmin, vmax)
im = ax.matshow(p0,
                norm = MidPointLogNorm(midpoint=1, vmin=vmin, vmax=vmax, log=np.log2),
                cmap = bwr,
                extent = extent)
plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
             ax=ax, ticks=locator2, format=formatter2)
ax.set_title("$p = 0$", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

# pseudocount = 0.0001
ax = axs[0, 1]
vmin = min(np.nanmin(p1), 1/np.nanmax(p1))
vmax = max(np.nanmax(p1), 1/np.nanmin(p1))
im = ax.matshow(p1,
                norm = MidPointLogNorm(midpoint=1, vmin=vmin, vmax=vmax, log=np.log2),
                cmap = bwr,
                extent = extent)
plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
             ax=ax, ticks=locator2, format=formatter2)
ax.set_title("$p = 0.0001$", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

# pseudocount = 0.001
ax = axs[1, 0]
vmin = min(np.nanmin(p2), 1/np.nanmax(p2))
vmax = max(np.nanmax(p2), 1/np.nanmin(p2))
im = ax.matshow(p2,
                norm = MidPointLogNorm(midpoint=1, vmin=vmin, vmax=vmax, log=np.log2),
                cmap = bwr,
                extent=extent)
plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction difference', 
             ax=ax, ticks=locator2, format=formatter2)
ax.set_title("$p = 0.001$", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

# pseudocount = 1
ax = axs[1, 1]
vmin = min(np.nanmin(p3), 1/np.nanmax(p3))
vmax = max(np.nanmax(p3), 1/np.nanmin(p3))
im = ax.matshow(p3,
                norm = MidPointLogNorm(midpoint=1, vmin=vmin, vmax=vmax, log=np.log2),
                cmap = bwr,
                extent = extent)
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, 
                    label='Interaction difference',
                    ax=ax, ticks=locator2, format=formatter2)
#cbar.ax.yaxis.set_major_formatter(mpl.ticker.NullFormatter())
cbar.ax.yaxis.set_major_formatter(mpl.ticker.LogFormatterMathtext(base=2, labelOnlyBase=True))
cbar.ax.yaxis.set_minor_formatter(mpl.ticker.LogFormatterMathtext(base=2, ))

ax.set_title("$p = 1$", y=1.01)

pplot.format_ticks(ax)
pplot.background_color(ax, 'gray-light')

plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', PATH_IMG, 
                         REGION + '_pseudocounts_log2_sym.pdf'))
