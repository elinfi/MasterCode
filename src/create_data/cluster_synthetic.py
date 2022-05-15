import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/notebooks/clustering')

import numpy as np

from cluster_functions import clustering_interaction


PATH_IN_DATA = '/home/elinfi/MasterCode/data/simulations/comparison'
PATH_OUT_DATA = '/home/elinfi/MasterCode/data/simulations/cluster/tad_2_3_0.5_tadtad_2_stripe_0.5'

REGION = 'chr10:6351511-10351511'
#EXTENSION = '_k_2'
#EXTENSION = '_tad_1.4_2_1.7'
#EXTENSION = '_tad_2_3_4_stripe_2'
EXTENSION = '_tad_2_3_0.5_tadtad_2_stripe_0.5'

MEDOIDS = [5, 6]
#RANDOM_STATE = [19, 45, 18888, 4321, 1791, 345628]
#RANDOM_STATE = [19]
RANDOM_STATE = [19, 922, 293, 473, 892, 113, 272, 249,  89, 458]

################################################################################
# CLUSTERING ###################################################################
################################################################################

for random_state in RANDOM_STATE:
    
    ############################################################################
    # SUBTRACTION ##############################################################
    ############################################################################

    sub = np.load(os.path.join(PATH_IN_DATA, 
                               REGION + EXTENSION + '_sub.npy'))
    clustering_interaction(sub, 
                           MEDOIDS,
                           PATH_OUT_DATA, 
                           REGION + EXTENSION + '_sub',
                           '',
                           random_state)

    ############################################################################
    # REALTIVE DIFFERENCE ######################################################
    ############################################################################

    reldiff = np.load(os.path.join(PATH_IN_DATA, 
                                   REGION + EXTENSION + '_reldiff.npy'))
    clustering_interaction(reldiff, 
                           MEDOIDS,
                           PATH_OUT_DATA, 
                           REGION + EXTENSION + '_reldiff',
                           '',
                           random_state)

    ############################################################################
    # POISSON ##################################################################
    ############################################################################

    poisson = np.load(os.path.join(PATH_IN_DATA, 
                                   REGION + EXTENSION + '_poisson.npy'))
    clustering_interaction(poisson, 
                           MEDOIDS,
                           PATH_OUT_DATA, 
                           REGION + EXTENSION + '_poisson',
                           '',
                           random_state)

    ############################################################################
    # LFC ######################################################################
    ############################################################################

    lfc = np.load(os.path.join(PATH_IN_DATA, 
                               REGION + EXTENSION + '_ratio.npy'))
    clustering_interaction(np.log2(lfc), 
                           MEDOIDS,
                           PATH_OUT_DATA, 
                           REGION + EXTENSION + '_lfc',
                           '',
                           random_state)

    ############################################################################
    # PSEUDO 0.0001, 0.001, 1 ##################################################
    ############################################################################

    p0 = np.load(os.path.join(PATH_IN_DATA, 
                               REGION + EXTENSION + f'_ratio.npy'))
    clustering_interaction(np.log2(p0), 
                           MEDOIDS,
                           PATH_OUT_DATA, 
                           REGION + EXTENSION + '_p0',
                           '',
                           random_state)

    for p in [0.0001, 0.001, 1]:
        pseudo = np.load(os.path.join(PATH_IN_DATA, 
                                      REGION + EXTENSION + f'_p_{p}.npy'))
        clustering_interaction(np.log2(pseudo), 
                               MEDOIDS,
                               PATH_OUT_DATA, 
                               REGION + EXTENSION + f'_p_{p}',
                               '',
                               random_state)