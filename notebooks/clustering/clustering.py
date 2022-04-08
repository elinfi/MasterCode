import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/class/')

import numpy as np

from cluster_functions import clustering


DATA_PATH = '/home/elinfi/MasterCode/data/simulations'
CLUSTER_PATH = '/home/elinfi/MasterCode/data/simulations/clusters'
REGION = 'chr10:6351511-10351511'
INTERACTION_WEIGHT = 1
DIAG_WEIGHT = [0, 0.001, 0.01, 0.1, 1]
MEDOIDS = 5
RANDOM_STATE = 19

# subtraction
sub = np.load(os.path.join(DATA_PATH, REGION + '_sub_k_2.npy'))
clustering(sub, INTERACTION_WEIGHT, DIAG_WEIGHT, MEDOIDS, CLUSTER_PATH, 
           f'{REGION}_sub', RANDOM_STATE)

# relative difference
rel = np.load(os.path.join(DATA_PATH, REGION + '_reldiff_k_2.npy'))
clustering(rel, INTERACTION_WEIGHT, DIAG_WEIGHT, MEDOIDS, CLUSTER_PATH, 
           f'{REGION}_reldiff', RANDOM_STATE)

# log fold change
ratio = np.load(os.path.join(DATA_PATH, REGION + '_ratio_k_2.npy'))
ratio_log = np.log2(ratio)
ratio_higlass = np.log2(ratio + np.nanmin(ratio[ratio > 0]))
clustering(ratio_log, INTERACTION_WEIGHT, DIAG_WEIGHT, MEDOIDS, CLUSTER_PATH, 
           f'{REGION}_ratio_log', RANDOM_STATE)
clustering(ratio_higlass, INTERACTION_WEIGHT, DIAG_WEIGHT, MEDOIDS, 
           CLUSTER_PATH, f'{REGION}_ratio_higlass', RANDOM_STATE)

# pseudocount
p001 = np.load(os.path.join(DATA_PATH, REGION + '_p001_k_2.npy'))
clustering(np.log10(p001), INTERACTION_WEIGHT, DIAG_WEIGHT, MEDOIDS, CLUSTER_PATH, 
           f'{REGION}_p001', RANDOM_STATE)
p01 = np.load(os.path.join(DATA_PATH, REGION + '_p01_k_2.npy'))
clustering(np.log10(p01), INTERACTION_WEIGHT, DIAG_WEIGHT, MEDOIDS, CLUSTER_PATH, 
           f'{REGION}_p01', RANDOM_STATE)
p1 = np.load(os.path.join(DATA_PATH, REGION + '_p1_k_2.npy'))
clustering(np.log10(p1), INTERACTION_WEIGHT, DIAG_WEIGHT, MEDOIDS, CLUSTER_PATH, 
           f'{REGION}_p1', RANDOM_STATE)

# poisson
poisson = np.load(os.path.join(DATA_PATH, REGION + '_poisson_k_2.npy'))
clustering(poisson, INTERACTION_WEIGHT, DIAG_WEIGHT, MEDOIDS, CLUSTER_PATH, 
           f'{REGION}_poisson', RANDOM_STATE)