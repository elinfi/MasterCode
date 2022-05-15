import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/class/')
sys.path.insert(2, '/home/elinfi/MasterCode/plot/func')
sys.path.insert(3, '/home/elinfi/MasterCode/notebooks/clustering')

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pretty_plotting as pplot

from rand_index import rand_index

PATH_IN = '/home/elinfi/MasterCode/data/simulations/cluster/tad_2_3_0.5_tadtad_2_stripe_0.5/'
PATH_OUT = '/home/elinfi/MasterCode/img/statistics/'
REGION = 'chr10:6351511-10351511'
SYNTHETIC = 'tad_2_3_0.5_tadtad_2_stripe_0.5'
EXTENSION = 'cluster_wdiag_0_medoids_4_rstate_19'
PATH_EXACT = '/home/elinfi/MasterCode/data/simulations/chr10:6351511-10351511_tad_2_3_0_tadtad_2_stripe_0_fasit.npy'

METHOD = ['sub', 'reldiff', 'lfc', 'p_0.0001', 'p_0.001', 'p_1']

NAMES = ['sub', 'reldiff', '$p = 0$', '$p = 10^{-4}$', '$p = 10^{-3}$', '$p = 1$', 'true']

################################################################################
# EXTRACT DATA #################################################################
################################################################################

list_filenames = []


# add estimated clusters to list_filenames
for method in METHOD:
    file = os.path.join(PATH_IN, f'{REGION}_{SYNTHETIC}_{method}_{EXTENSION}.npy')
    list_filenames.append(file)
    
# add true clusters to list_filenames
list_filenames.append(PATH_EXACT)
    
# get rand index score
rand_idxes = rand_index(list_filenames)

################################################################################
# STYLE SETTINGS ###############################################################
################################################################################

# use seaborn style
sns.set_theme('paper')
sns.set_style('ticks')

# set font sizes
pplot.font_size(16, 18, 20)

################################################################################
# PLOT RAND INDEX ##############################################################
################################################################################

f, axs = plt.subplots(figsize=(8, 6.5),
                      nrows=1,
                      ncols=1,
                      sharex=True, sharey=False)

# adjust subplots
plt.subplots_adjust(top=0.95,
                    bottom=0.15,
                    left=0.14,
                    right=0.999)

ax = axs

final_rand_idxes = rand_idxes[1:, :-1]
mask = np.zeros_like(final_rand_idxes)
mask[np.triu_indices_from(mask, 1)] = True

sns.heatmap(final_rand_idxes, 
            vmin=0, vmax=1, 
            linewidths=.1, 
            annot=True, square=False, mask=mask,
            xticklabels=NAMES[:-1],
            yticklabels=NAMES[1:],
            ax=ax,
            cbar_kws={'pad': 0.02,
                      'label': 'rand index score'})
ax.set_title('Rand index', y=0)
plt.xticks(rotation=35)
plt.yticks(rotation=0)

plt.savefig(os.path.join(PATH_OUT, f'{REGION}_{SYNTHETIC}_{EXTENSION}_randidx_adjusted.pdf'))
