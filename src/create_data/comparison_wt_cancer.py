import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/func')

import numpy as np
import comparison as comp


REGION = 'chr12:115900000-120200000'
RESOLUTION = 16000
PATH_IN_DATA = '/home/elinfi/MasterCode/data/wt_cancer/raw/'
PATH_OUT_DATA = '/home/elinfi/MasterCode/data/wt_cancer/comparison'

PATH_WT = os.path.join(PATH_IN_DATA, 
                       f'{REGION}_res_{RESOLUTION}_wt_merged.npy')
PATH_CANCER = os.path.join(PATH_IN_DATA, 
                           f'{REGION}_res_{RESOLUTION}_cancer_merged.npy')


################################################################################
# PSEUDOCOUNT 0.0001, 0.001, 1 #################################################
################################################################################
wt = np.load(PATH_WT)
cancer = np.load(PATH_CANCER)

for pseudo in [0.0001, 0.001, 1]:
    p = comp.pseudo_diff(wt, cancer, pseudo)
    np.save(os.path.join(PATH_OUT_DATA, 
                         f'{REGION}_res_{RESOLUTION}_wt_cancer_p_{pseudo}.npy'),
            p)

################################################################################
# SUBTRACTION ##################################################################
################################################################################
wt = np.load(PATH_WT)
cancer = np.load(PATH_CANCER)

sub = comp.subtraction(wt, cancer)
np.save(os.path.join(PATH_OUT_DATA, 
                     f'{REGION}_res_{RESOLUTION}_wt_cancer_sub.npy'),
        sub)

################################################################################
# RELDIFF ######################################################################
################################################################################
wt = np.load(PATH_WT)
cancer = np.load(PATH_CANCER)

reldiff = comp.relative_difference(wt, cancer)
np.save(os.path.join(PATH_OUT_DATA, 
                     f'{REGION}_res_{RESOLUTION}_wt_cancer_reldiff.npy'),
        reldiff)

################################################################################
# POISSON ######################################################################
################################################################################
wt = np.load(PATH_WT)
cancer = np.load(PATH_CANCER)

poisson = comp.poisson(wt, cancer)
np.save(os.path.join(PATH_OUT_DATA, 
                     f'{REGION}_res_{RESOLUTION}_wt_cancer_poisson.npy'),
        poisson)

################################################################################
# FOLD CHANGE ##################################################################
################################################################################
wt = np.load(PATH_WT)
cancer = np.load(PATH_CANCER)

ratio = comp.fold_change(wt, cancer)
np.save(os.path.join(PATH_OUT_DATA, 
                     f'{REGION}_res_{RESOLUTION}_wt_cancer_ratio.npy'),
        ratio)
