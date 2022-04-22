import os

import numpy as np

from k_medoids import KMedoids
from dissimilarity_matrix import DissimilarityMatrix

def clustering(data, interaction_weight, diag_weight, medoids, path, filename, 
               random_state):
    diss_data = DissimilarityMatrix(data)
    n = diss_data.n
    triu_nan_idx = diss_data.nan_idx_triu
    
    # interactions difference
    interaction_dist = diss_data.scipy_dist(metric='interactions_dist', 
                                           col1=0, col2=3)
    
    # diagonal distance
    diag_dist = diss_data.scipy_dist(metric='diagonal_dist', col1=0, col2=3)
    
    if type(diag_weight) == list:
        for w2 in diag_weight:
            dist_mat = interaction_weight*interaction_dist + w2*diag_dist
            
            kmedoids = KMedoids(dist_mat)
            cluster_result = kmedoids.clusters(medoids=medoids, 
                                              random_state=random_state)
            labels_mat = kmedoids.labels4plotting_nan(n, triu_nan_idx)
            
            np.save(os.path.join(path, 
                                 filename 
                                 + f'_cluster_wdiag_{w2}_medoids_{medoids}.npy'), 
                    labels_mat)
    else:
        dist_mat = interaction_weight*interaction_dist + diag_weight*diag_dist
            
        kmedoids = KMedoids(dist_mat)
        cluster_result = kmedoids.clusters(medoids=medoids, 
                                           random_state=random_state)
        labels_mat = kmedoids.labels4plotting_nan(n, triu_nan_idx)

        np.save(os.path.join(path, 
                             filename 
                             + f'_cluster_wdiag_{diag_weight}_medoids_{medoids}.npy'), 
                labels_mat)
        
def clustering_interaction(data, medoids, path, filename, extension, random_state):
    
    diss_data = DissimilarityMatrix(data)
    n = diss_data.n
    triu_nan_idx = diss_data.nan_idx_triu
    
    # interactions difference
    dist_mat = diss_data.scipy_dist(metric='interactions_dist', col1=0, col2=3)

    kmedoids = KMedoids(dist_mat)
    cluster_result = kmedoids.clusters(medoids=medoids, 
                                       random_state=random_state)
    labels_mat = kmedoids.labels4plotting_nan(n, triu_nan_idx)

    np.save(os.path.join(path, 
                         filename 
                         + f'_cluster_wdiag_0_medoids_{medoids}{extension}_r_{random_state}.npy'), 
            labels_mat)
        
        
def save_distance_matrices(data, path, filename):
    diss_data = DissimilarityMatrix(data)
    n = diss_data.n
    triu_nan_idx = diss_data.nan_idx_triu
    
    # interactions difference
    interaction_dist = diss_data.scipy_dist(metric='interactions_dist', 
                                           col1=0, col2=3)
    
    # diagonal distance
    diag_dist = diss_data.scipy_dist(metric='diagonal_dist', col1=0, col2=3)
    
    with open(os.join.path(path, filename + '_dist_mat.npy')) as f:
              np.save(f, n)
              np.save(f, triu_nan_idx)
              np.save(f, interaction_dist)
              np.save(f, diag_dist)
              
def open_distance_matrices(distmat_file, interaction_weight, diag_weight, 
                           medoids, path, filename, random_state):
    """
    Open file with distance matrices and performs clustering.
    
    Args:
        distmat_file (string):
            path to file containing distance matrices
        interaction_weight (float):
            weight for pixel value (should always be 1)
        diag_weight (float/list):
            weight for diagonal distance. float or list of floats.
        medoids (int):
            number of medoids.
        path (string):
            path to where to store data
        filename (string):
            beginning of filename on the form '{region}'
        random_state (int):
            random state 
    """
              
    with open(distmat_file) as f:
              n = np.load(f)
              triu_nan_idx = np.load(f)
              interaction_dist = np.load(f)
              diag_dist = np.load(f)
              
    if type(diag_weight) == list:
        for w2 in diag_weight:
            dist_mat = interaction_weight*interaction_dist + w2*diag_dist
            
            kmedoids = KMedoids(dist_mat)
            cluster_result = kmedoids.clusters(medoids=medoids, 
                                              random_state=random_state)
            labels_mat = kmedoids.labels4plotting_nan(n, triu_nan_idx)
            
            np.save(os.path.join(path, filename + f'_cluster_wdiag_{w2}_medoids_{medoids}.npy'), 
                    labels_mat)
    else:
        dist_mat = interaction_weight*interaction_dist + diag_weight*diag_dist
            
        kmedoids = KMedoids(dist_mat)
        cluster_result = kmedoids.clusters(medoids=medoids, 
                                          random_state=random_state)
        labels_mat = kmedoids.labels4plotting_nan(n, triu_nan_idx)

        np.save(os.path.join(path, 
                             filename 
                             + f'_cluster_wdiag_{diag_weight}_medoids_{medoids}.npy'), 
                labels_mat)