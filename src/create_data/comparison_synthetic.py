import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/func')

import numpy as np
import comparison as comp


REGION = 'chr10:6351511-10351511'
#EXTENSION = '_tad_1.4_2_1.7'
#EXTENSION = '_k_2'
EXTENSION = '_tad_2_3_4_tadtad_3_stripe_2'
PATH_IN_DATA = '/home/elinfi/MasterCode/data/simulations/'
PATH_WT1 = os.path.join(PATH_IN_DATA, REGION + '_wt1.npy')
PATH_MOD = os.path.join(PATH_IN_DATA, REGION + EXTENSION + '.npy')
PATH_OUT_DATA = '/home/elinfi/MasterCode/data/simulations/comparison'

################################################################################
# LOAD SYNTHETIC DATA ##########################################################
################################################################################

wt1 = np.load(PATH_WT1)
mod = np.load(PATH_MOD)

################################################################################
# PSEUDOCOUNT 0.0001, 0.001, 1 #################################################
################################################################################
wt1 = np.load(PATH_WT1)
mod = np.load(PATH_MOD)

p1 = comp.pseudo_diff(wt1, mod, 0.0001)
p2 = comp.pseudo_diff(wt1, mod, 0.001)
p3 = comp.pseudo_diff(wt1, mod, 1)

print(np.sum(np.isnan(p1)))
print(np.sum(np.isnan(p2)))
print(np.sum(np.isnan(p3)))

np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + f'_p_{0.0001}.npy'),
        p1)
np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + f'_p_{0.001}.npy'),
        p2)
np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + f'_p_{1}.npy'),
        p3)

################################################################################
# SUBTRACTION ##################################################################
################################################################################
wt1 = np.load(PATH_WT1)
mod = np.load(PATH_MOD)

sub = comp.subtraction(wt1, mod)
np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + '_sub.npy'),
        sub)

################################################################################
# RELDIFF ######################################################################
################################################################################
wt1 = np.load(PATH_WT1)
mod = np.load(PATH_MOD)

reldiff = comp.relative_difference(wt1, mod)
np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + '_reldiff.npy'),
        reldiff)

################################################################################
# POISSON ######################################################################
################################################################################
wt1 = np.load(PATH_WT1)
mod = np.load(PATH_MOD)

poisson = comp.poisson(wt1, mod)
np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + '_poisson.npy'),
        poisson)

################################################################################
# FOLD CHANGE ##################################################################
################################################################################
wt1 = np.load(PATH_WT1)
mod = np.load(PATH_MOD)

ratio = comp.fold_change(wt1, mod)
np.save(os.path.join(PATH_OUT_DATA, REGION + EXTENSION + '_ratio.npy'),
        ratio)
