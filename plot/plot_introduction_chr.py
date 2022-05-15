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
from data_preparation import DataPreparation

REGION = 'chr4:0-40000000'
RESOLUTION = 16000
BALANCE = True
PATH_WT = '/home/elinfi/coolers/HiC_wt_merged.mcool'
PATH_IMG = '/home/elinfi/MasterCode/img/'

# create objects of class
wt = DataPreparation(PATH_WT, RESOLUTION, REGION, BALANCE)


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

# create subplots
f, axs = plt.subplots(figsize=(9.1, 7.7),
                      nrows=1,
                      ncols=1,
                      sharex=True, sharey=True)

# adjust subplots
plt.subplots_adjust(left=0.01,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.95)

# wt 001
ax = axs
im = ax.matshow(wt.matrix,
                cmap=fall,
                norm= LogNorm(),
                extent=extent
               )
ax.set_title('Chromosome 10', y=1.01)
#plt.setp(ax.get_xticklabels(), visible=False)
#plt.setp(ax.get_yticklabels(), visible=False)
#ax.tick_params(axis='both', which='both', length=0)

ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
pplot.format_ticks(ax)


plt.colorbar(im, fraction=0.046, pad=0.04, label='Interaction frequency',
             ticks=LogLocator(base=2), 
             format=LogFormatterSciNotation(base=2))
    
plt.savefig(os.path.join(PATH_IMG, 
                         f'{REGION}_res_{RESOLUTION}_wt_merged.pdf'))
