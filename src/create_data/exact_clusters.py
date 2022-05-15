import sys
sys.path.append('/home/elinfi/MasterCode/src/class/')
sys.path.append('/home/elinfi/MasterCode/src/func/')
sys.path.append('/home/elinfi/MasterCode/notebooks/plotting/')
sys.path.append('/home/elinfi/MasterCode/plot/func')

import cooltools.lib.plotting

import numpy as np
import matplotlib as mpl
import pretty_plotting as pplot
import plot_subtraction as sub
import plot_differences as splot
import matplotlib.pyplot as plt
import change_functions as cf

from exact_clustering import ExactClustering
from mid_point_log_norm import MidPointLogNorm
from fplot_replicates import plot_replicates

data = ExactClustering()
data.print_tads_loops()
print(np.sum((data.org + data.mod) == data.org)/(data.org.shape[0]*data.org.shape[1]))
print(np.sum(np.isnan(data.mod)))

k = 1
k1 = 2    # tad
k2 = 3    # tad
k3 = 1/2  # tad
k4 = k    # dot
k5 = k    # dot
k6 = 2    # tad-tad
k7 = 1/2    # stripe


data.change_tad(cf.constant, [0], k=k1)
data.change_tad(cf.constant, [2], k=k2)
data.change_tad(cf.constant, [3], k=k3)
data.change_loop(cf.constant, [0], k=k4)
data.change_loop(cf.constant, [1], k=k5)
data.change_tad_tad(cf.constant, [0, 2], k=k6)
data.change_stripe(5, 157, 186, cf.constant, k=k7)

np.save(f'/home/elinfi/MasterCode/data/simulations/{data.mat_region}_tad_{k1}_{k2}_{k3}_tadtad_{k6}_stripe_{k7}_exact.npy', data.mod)
#np.save(f'/home/elinfi/MasterCode/data/simulations/{data.mat_region}_wt1.npy', data.org)