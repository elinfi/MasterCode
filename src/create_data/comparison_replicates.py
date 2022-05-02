import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/func')

import numpy as np
import comparison as comp


REGION = 'chr10:6351511-10351511'
EXTENSION = ''
PATH_IN_DATA = '/home/elinfi/MasterCode/data/replicates/'
PATH_WT1 = os.path.join(PATH_IN_DATA, REGION + '_wt1.npy')
PATH_WT2 = os.path.join(PATH_IN_DATA, REGION + '_wt2.npy')
PATH_OUT_DATA = '/home/elinfi/MasterCode/data/replicates/comparison'

################################################################################
# PSEUDOCOUNT 0.0001, 0.001, 1 #################################################
################################################################################
wt1 = np.load(PATH_WT1)
wt2 = np.load(PATH_WT2)

for pseudo in [0.0001, 0.001, 0.01, 0.1, 1]:
    p = comp.pseudo_diff(wt1, wt2, pseudo)
    np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + f'_p_{pseudo}.npy'),
            p)

################################################################################
# SUBTRACTION ##################################################################
################################################################################
wt1 = np.load(PATH_WT1)
wt2 = np.load(PATH_WT2)

sub = comp.subtraction(wt1, wt2)
np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + '_sub.npy'),
        sub)

################################################################################
# RELDIFF ######################################################################
################################################################################
wt1 = np.load(PATH_WT1)
wt2 = np.load(PATH_WT2)

reldiff = comp.relative_difference(wt1, wt2)
np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + '_reldiff.npy'),
        reldiff)

################################################################################
# POISSON ######################################################################
################################################################################
wt1 = np.load(PATH_WT1)
wt2 = np.load(PATH_WT2)

poisson = comp.poisson(wt1, wt2)
np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + '_poisson.npy'),
        poisson)

################################################################################
# FOLD CHANGE ##################################################################
################################################################################
wt1 = np.load(PATH_WT1)
wt2 = np.load(PATH_WT2)

ratio = comp.fold_change(wt1, wt2)
np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + '_ratio.npy'),
        ratio)
