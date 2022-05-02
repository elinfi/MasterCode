"""
Saves dissimilarity matrix for synthetic data. The datafiles produced are too
large, thus this method should not be used.
"""

import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/notebooks/clustering')

import numpy as np

from cluster_functions import save_distance_matrices

REGION = 'chr10:6351511-10351511'
#EXTENSION = '_tad_1.2_2_1.7'
EXTENSION = '_k_2'
PATH_IN_DATA = '/home/elinfi/MasterCode/data/simulations/comparison'
PATH_OUT_DATA = '/home/elinfi/MasterCode/data/simulations/dissimilarity'

################################################################################
# SUBTRACTION ##################################################################
################################################################################

sub = np.load(os.path.join(PATH_IN_DATA, REGION + EXTENSION + '_sub.npy'))
save_distance_matrices(sub, PATH_OUT_DATA, REGION + EXTENSION + '_sub')

################################################################################
# REALTIVE DIFFERENCE ##########################################################
################################################################################

reldiff = np.load(os.path.join(PATH_IN_DATA, REGION + EXTENSION + '_reldiff.npy'))
save_distance_matrices(reldiff, PATH_OUT_DATA, REGION + EXTENSION + '_reldiff')

################################################################################
# POISSON ######################################################################
################################################################################

poisson = np.load(os.path.join(PATH_IN_DATA, REGION + EXTENSION + '_poisson.npy'))
save_distance_matrices(poisson, PATH_OUT_DATA, REGION + EXTENSION + '_poisson')

################################################################################
# LFC ##########################################################################
################################################################################

lfc = np.load(os.path.join(PATH_IN_DATA, REGION + EXTENSION + '_ratio.npy'))
save_distance_matrices(np.log2(lfc), PATH_OUT_DATA, REGION + EXTENSION + '_lfc')

################################################################################
# PSEUDO 0.0001, 0.001, 1 #####################################################
################################################################################

p0 = np.load(os.path.joing(PATH_IN_DATA, REGION + EXTENSION + f'_ratio.npy'))
save_distance_matrices(np.log10(p0), PATH_OUT_DATA, REGION + EXTENSION + '_p_0')

for p in [0.0001, 0.001, 1]:
    pseudo = np.load(os.path.join(PATH_IN_DATA, 
                                  REGION + EXTENSION + f'_p_{p}.npy'))
    save_distance_matrices(np.log10(pseudo), 
                         PATH_OUT_DATA, REGION + EXTENSION + f'_p_{p}')





