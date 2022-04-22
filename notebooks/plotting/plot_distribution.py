import os
import sys
import cooltools.lib.plotting
sys.path.insert(1, '/home/elinfi/MasterCode/src/class/')

import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pretty_plotting as pplot

from copy import copy
from matplotlib.ticker import LogLocator, LogFormatterSciNotation
from mid_point_log_norm import MidPointLogNorm
from matplotlib.colors import TwoSlopeNorm, LogNorm
from sklearn.preprocessing import StandardScaler, normalize

PATH_DIFF = '/home/elinfi/MasterCode/data/replicates/'
PATH_IMG = 'replicates/pdf/distribution'
REGION = 'chr10:6351511-10351511'

################################################################################
# STYLE SETTINGS ###############################################################
################################################################################

# use seaborn style
sns.set_theme('paper')
sns.set_style('ticks')

# set font sizes
pplot.font_size(16, 18, 20)

################################################################################
# PSEUDOCOUNT 0, 0.0001, 0.001, 1 ##############################################
################################################################################

p0 = np.load(os.path.join(PATH_DIFF, f'{REGION}_ratio.npy'))
p1 = np.load(os.path.join(PATH_DIFF, f'{REGION}_p00001.npy'))
p2 = np.load(os.path.join(PATH_DIFF, f'{REGION}_p0001.npy'))
p3 = np.load(os.path.join(PATH_DIFF, f'{REGION}_p1.npy'))

# only use upper triangular data because of symmetry
n = p0.shape[0]
triu = np.triu_indices(n)

p0 = p0[triu]
p1 = p1[triu]
p2 = p2[triu]
p3 = p3[triu]

p0 = p0[p0 > 0]
"""
p0 = p0[~np.isnan(p0)]
p1 = p1[~np.isnan(p1)]
p2 = p2[~np.isnan(p2)]
p3 = p3[~np.isnan(p3)]


# normalize
p0 = normalize(p0.reshape(1, -1))
p1 = normalize(p1.reshape(1, -1))
p2 = normalize(p2.reshape(1, -1))
p3 = normalize(p3.reshape(1, -1))
"""
# create subplot
fig, axs = plt.subplots(figsize=(16, 12),
                        nrows=2,
                        ncols=2,
                        sharex=False, sharey=False)
# adjust subplots
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.94,
                    wspace=0.3)

sns.histplot(data=p0, ax=axs[0, 0], bins=50, log_scale=True)
sns.histplot(data=p1, ax=axs[0, 0], bins=50, log_scale=True)
sns.histplot(data=p2, ax=axs[1, 0], bins=50, log_scale=True)
sns.histplot(data=p3, ax=axs[1, 1], bins=50, log_scale=True)
plt.setp(axs[1,1].get_xticklabels(), rotation=45, ha='right')
"""
sns.kdeplot(data=p0, ax=axs[0, 0])
sns.kdeplot(data=p1, ax=axs[0, 1])
sns.kdeplot(data=p2, ax=axs[1, 0])
sns.kdeplot(data=p3, ax=axs[1, 1])
"""
#penguins = sns.load_dataset("penguins")
#sns.displot(data=penguins, x="flipper_length_mm")
#print(penguins)
plt.savefig(os.path.join('/home/elinfi/MasterCode/img/', PATH_IMG, 
                         f'{REGION}_distribution_pseudocount_log.pdf'))