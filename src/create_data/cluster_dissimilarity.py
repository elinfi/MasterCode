import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/notebooks/clustering')

import numpy as np

REGION = 'chr10:6351511-10351511'
#EXTENSION = '_tad_1.2_2_1.7'
EXTENSION = '_k_2'
PATH_IN_DATA = '/home/elinfi/MasterCode/data/simulations/dissimilarity'
PATH_OUT_DATA = '/home/elinfi/MasterCode/data/simulations/clusters'

MEDOIDS = 2
INTERACTION_WEIGHT = 1
DIAG_WEIGHT = 0
RANDOM_STATE = 19

################################################################################
# SUBTRACTION ##################################################################
################################################################################

sub_file = os.path.join(PATH_IN_DATA, REGION + EXTENSION + '_sub_dist_mat.npy')
open_distance_matrices(sub_file, PATH_OUT_DATA, REGION + EXTENSION + '_sub')

################################################################################
# REALTIVE DIFFERENCE ##########################################################
################################################################################

reldiff = np.load(os.path.join(PATH_IN_DATA, 
                               REGION + EXTENSION + '_reldiff_dist_mat.npy'))
open_distance_matrices(reldiff, 
                       INTERACTION_WEIGHT,
                       DIAG_WEIGHT,
                       PATH_OUT_DATA, 
                       REGION + EXTENSION + '_reldiff',
                       RANDOM_STATE)

################################################################################
# POISSON ######################################################################
################################################################################

poisson = np.load(os.path.join(PATH_IN_DATA, 
                               REGION + EXTENSION + '_poisson_dist_mat.npy'))
open_distance_matrices(poisson, PATH_OUT_DATA, REGION + EXTENSION + '_poisson')

################################################################################
# LFC ##########################################################################
################################################################################

lfc = np.load(os.path.join(PATH_IN_DATA, 
                           REGION + EXTENSION + '_ratio_dist_mat.npy'))
open_distance_matrices(np.log2(lfc), PATH_OUT_DATA, REGION + EXTENSION + '_lfc')

################################################################################
# PSEUDO 0.0001, 0.001, 1 ######################################################
################################################################################

p0 = np.load(os.path.joing(PATH_IN_DATA, 
                           REGION + EXTENSION + f'_p_0_dist_mat.npy'))
open_distance_matrices(np.log10(p0), PATH_OUT_DATA, REGION + EXTENSION + '_p_0')

for p in [0.0001, 0.001, 1]:
    pseudo = np.load(os.path.join(PATH_IN_DATA, 
                                  REGION + EXTENSION + f'_p_{p}_dist_mat.npy'))
    open_distance_matrices(np.log10(pseudo), 
                         PATH_OUT_DATA, REGION + EXTENSION + f'_p_{p}')





