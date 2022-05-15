import os

import numpy as np

PATH = '/home/elinfi/MasterCode/data/simulations/cluster/tad_2_3_0.5_tadtad_2_stripe_0.5/'
REGION = 'chr10:6351511-10351511'
SYNTHETIC = 'tad_2_3_0.5_tadtad_2_stripe_0.5'
EXTENSION = 'cluster_wdiag_0_medoids_4'

METHOD = ['sub', 'reldiff', 'poisson', 'lfc', 'p0', 'p_0.0001', 'p_0.001', 'p_1']
RSTATE = [19, 922, 293, 473, 892, 113, 272, 249,  89, 458]

for method in METHOD:
    matrix = []
    for rstate in RSTATE:
        filename = os.path.join(PATH, f'{REGION}_{SYNTHETIC}_{method}_{EXTENSION}_rstate_{rstate}.npy')
        matrix.append(np.load(filename))
    print('hei')
    print(np.unique(np.array(matrix)))
    median = np.round(np.nanmean(np.array(matrix), axis=0))
    print(np.unique(median))
    file_out = os.path.join(PATH, f'{REGION}_{SYNTHETIC}_{method}_{EXTENSION}_median.npy')
    np.save(file_out, median)