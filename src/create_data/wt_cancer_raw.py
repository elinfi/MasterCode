import os
import sys
sys.path.insert(1, '/home/elinfi/MasterCode/src/class')

import numpy as np

from data_preparation import DataPreparation

REGION = 'chr8:108000000-112000000'
RESOLUTION = 16000
BALANCE = True
EXTENSION = ''
PATH_WT = '/home/elinfi/coolers/HiC_wt_merged.mcool'
PATH_CANCER = '/home/elinfi/coolers/HiC_cancer_merged.mcool'
PATH_OUT_DATA = '/home/elinfi/MasterCode/data/wt_cancer/'

# create objects of class
wt = DataPreparation(PATH_WT, RESOLUTION, REGION, BALANCE)
cancer = DataPreparation(PATH_CANCER, RESOLUTION, REGION, BALANCE)

np.save(os.path.join(PATH_OUT_DATA, REGION + '_wt_merged.npy'), wt.matrix)
np.save(os.path.join(PATH_OUT_DATA, REGION + '_cancer_merged.npy'), cancer.matrix)
import os
import sys
sys.path.append('/home/elinfi/MasterCode/src/class/')

from data_preparation import DataPreparation

PATH_WT = '/home/elinfi/coolers/HiC_wt_merged.mcool'
PATH_CANCER = '/home/elinfi/coolers/HiC_cancer_merged.mcool'
REGION = 'chr12:115900000-120200000'
RESOLUTION = 16000
BALANCE = True

PATH_OUT = '/home/elinfi/MasterCode/data/wt_cancer/raw/'


wt = DataPreparation(PATH_WT, RESOLUTION, REGION, BALANCE)
cancer = DataPreparation(PATH_CANCER, RESOLUTION, REGION, BALANCE)

# save data
np.save(os.path.join(PATH_OUT, f'{REGION}_res_{RESOLUTION}_wt_merged.npy'), 
        wt.matrix)
np.save(os.path.join(PATH_OUT, f'{REGION}_res_{RESOLUTION}_cancer_merged.npy'), 
        cancer.matrix)
