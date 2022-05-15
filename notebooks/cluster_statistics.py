import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/class/')
sys.path.insert(2, '/home/elinfi/MasterCode/plot/func')
sys.path.insert(3, '/home/elinfi/MasterCode/notebooks/clustering')

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pretty_plotting as pplot

from rand_index import rand_index

PATH_IN = '/home/elinfi/MasterCode/data/simulations/cluster/tad_2_3_0.5_tadtad_2_stripe_0.5/'
PATH_OUT = '/home/elinfi/MasterCode/img/statistics/'
REGION = 'chr10:6351511-10351511'
SYNTHETIC = 'tad_2_3_0.5_tadtad_2_stripe_0.5'
EXTENSION = 'cluster_wdiag_0_medoids_4_rstate_19'
PATH_EXACT = '/home/elinfi/MasterCode/data/simulations/chr10:6351511-10351511_tad_2_3_4_tadtad_5_stripe_6_fasit.npy'

true_clusters = np.load(PATH_EXACT)
tad1 = true_clusters == 2
tad2 = true_clusters == 3
tad3 = true_clusters == 4
tadtad = true_clusters == 5
stripe = true_clusters == 6

