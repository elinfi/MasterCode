import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/class/')

import numpy as np

from cluster_functions import clustering_interaction


DATA_PATH = '/home/elinfi/MasterCode/data/simulations'
CLUSTER_PATH = '/home/elinfi/MasterCode/data/simulations/clusters/k2'
EXTENSION = '_k_2'
REGION = 'chr10:6351511-10351511'
MEDOIDS = 2
RANDOM_STATE = [19, 45, 18888, 4321, 1791, 345628]

for i in range(len(RANDOM_STATE)):
    # subtraction
    sub = np.load(os.path.join(DATA_PATH, REGION + '_sub_k_2.npy'))
    clustering_interaction(sub, MEDOIDS, CLUSTER_PATH, 
               f'{REGION}_sub', EXTENSION, RANDOM_STATE[i])
    
    # relative difference
    rel = np.load(os.path.join(DATA_PATH, REGION + '_reldiff_k_2.npy'))
    clustering_interaction(rel, MEDOIDS, CLUSTER_PATH, 
               f'{REGION}_reldiff', EXTENSION, RANDOM_STATE[i])

    # log fold change
    ratio = np.load(os.path.join(DATA_PATH, REGION + '_ratio_k_2.npy'))
    ratio_log = np.log2(ratio)
    ratio_higlass = np.log2(ratio + np.nanmin(ratio[ratio > 0]))
    clustering_interaction(ratio_log, MEDOIDS, CLUSTER_PATH, 
               f'{REGION}_ratio_log', EXTENSION, RANDOM_STATE[i])
    clustering_interaction(ratio_higlass, MEDOIDS, 
               CLUSTER_PATH, f'{REGION}_ratio_higlass', EXTENSION, RANDOM_STATE[i])

    # pseudocount
    p001 = np.load(os.path.join(DATA_PATH, REGION + '_p001_k_2.npy'))
    clustering_interaction(np.log10(p001), MEDOIDS, CLUSTER_PATH, 
               f'{REGION}_p001', EXTENSION, RANDOM_STATE[i])
    p01 = np.load(os.path.join(DATA_PATH, REGION + '_p01_k_2.npy'))
    clustering_interaction(np.log10(p01), MEDOIDS, CLUSTER_PATH, 
               f'{REGION}_p01', EXTENSION, RANDOM_STATE[i])
    p1 = np.load(os.path.join(DATA_PATH, REGION + '_p1_k_2.npy'))
    clustering_interaction(np.log10(p1), MEDOIDS, CLUSTER_PATH, 
               f'{REGION}_p1', EXTENSION, RANDOM_STATE[i])

    # poisson
    poisson = np.load(os.path.join(DATA_PATH, REGION + '_poisson_k_2.npy'))
    clustering_interaction(poisson, MEDOIDS, CLUSTER_PATH, 
               f'{REGION}_poisson', EXTENSION, RANDOM_STATE[i])