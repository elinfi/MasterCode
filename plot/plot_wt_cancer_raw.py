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

REGION = 'chr12:115900000-120200000'
RESOLUTION = 16000
PATH_IN_DATA = '/home/elinfi/MasterCode/data/wt_cancer/raw/'
PATH_IMG = '/home/elinfi/MasterCode/img/wt_cancer/pdf/raw'

PATH_WT = os.path.join(PATH_IN_DATA, 
                       f'{REGION}_res_{RESOLUTION}_wt_merged.npy')
PATH_CANCER = os.path.join(PATH_IN_DATA, 
                           f'{REGION}_res_{RESOLUTION}_cancer_merged.npy')

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
# WT & CANCER ##################################################################
################################################################################

wt = np.load(PATH_WT)
cancer = np.load(PATH_CANCER)

# create subplots
f, axs = plt.subplots(figsize=(15.5, 7.7),
                      nrows=1,
                      ncols=2,
                      sharex=True, sharey=True)

# adjust subplots
plt.subplots_adjust(left=0.075,
                    bottom=0.08, 
                    right=0.94, 
                    top=0.99, 
                    wspace=0.04)

# wt 001
ax = axs[0]
im = ax.matshow(wt,
                cmap=fall,
                norm= LogNorm(),
                extent=extent)
ax.set_title('Wild type', y=1.01)
ax.set_title('A', loc='left', fontweight="bold", y=1.01)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
pplot.format_ticks(ax)

# wt 002
ax = axs[1]
im = ax.matshow(cancer,
                cmap=fall,
                norm=LogNorm(),
                extent=extent)

ax.set_title('Breast cancer', y=1.01)
ax.set_title('B', loc='left', fontweight="bold", y=1.01)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
pplot.format_ticks(ax)

plt.colorbar(im, fraction=0.0235, pad=0.02, label='Interaction frequency', 
             ax=list(axs), ticks=LogLocator(base=2), 
             format=LogFormatterSciNotation(base=2))
    
plt.savefig(os.path.join(PATH_IMG, 
                         f'{REGION}_res_{RESOLUTION}_wt_cancer.pdf'))