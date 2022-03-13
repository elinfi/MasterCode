
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/clustering/class/')

import numpy as np
import simulation_functions as sim

def test_tad_indexing():
    mat_region = 'chr8:85013332-95013332'
    tad_region = 'chr8:87000000-87600000'
    resolution = 32000
    replicate = 'wt_001'
    
    tad_mat = sim.matrix(tad_region, replicate, resolution)
    mat_mat = sim.matrix(mat_region, replicate, resolution)
    
    tad_i, tad_j = sim.get_tad_idx(mat_region, tad_region, resolution)
    
    tad_i = tad_i
    tad_j = tad_j
    
    tad_idx = mat_mat[tad_i:tad_j, tad_i:tad_j]
    
    test_bool = tad_mat != tad_idx
    test = np.sum(test_bool)
    print(tad_mat != tad_idx)
    
    print(tad_mat[test_bool])
    print(tad_idx[test_bool])
    
    if test == 0:
        return 'Hurray!'
    else:
        return f'It is wrong... Test = {test}'
    
print(test_tad_indexing())