import kmedoids

import numpy as np
import numpy.ma as ma

from numba import jit
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import MinMaxScaler



class DissimilarityMatrix:
    def __init__(self, data, remove_nan=True, scaler=None):
        """Creates a dissimilarity matrix.

        Creates a dissimilarity based on a Hi-C data matrix with a given
        distance metric.

        Args:
            data (ndarray):
                Square symmetric matrix containing the Hi-C data to cluster
            metric {'euclidean', 'test_euclidean'}:
                Which metric to use when calculating the dissimilarity matrix
            remove_nan (boolean):
                True when removing NaN's, False when keeping NaN's
            scaler ():
                
        Attributes:
            n (int):
                Length of data matrix
            triu_idx (array):
                indices for upper triangular matrix of n x n matrix
            nan_idx (boolean array):
                boolean matrix where all NaN's in data are true
            nan_idx_triu (array):
                upper triangular matrix of nan_idx
            X (ndarray):
                n x 3 array containing the bin1 and bin2 in the first two
                columns and the number of interactions in the last column
        """
        self.data = data
        self.remove_nan = remove_nan
        self.scaler = scaler
        
        # get shape of data matrix
        self.n = data.shape[0]
        
        # get indexes for upper triangular matrix of n x n matrix
        self.triu_idx = np.triu_indices(self.n)
        
        # get boolean matrix with true for all NaN's
        self.nan_idx = np.isnan(self.data)
        
        # get upper triangluar matrix of boolean NaN matrix
        self.nan_idx_triu = self.nan_idx[self.triu_idx]
        
        # transform data
        self.X = self.transform_data()
            
    def transform_data(self):
        """Extract the upper triangular of the data matrix and removes all NaN.

        Return:
            X (ndarray):
                n x 3 array containing the row indices in the first column,
                column indices in the second colum and the number of 
                interactions in the last column.
        """
        #assert log_scale in [None, 'log2', 'log10', 'MinMax'], 'Only None, 2 or 10 are valid log scales'
            
        n = self.n
        
        # create 2d coordinate system in 1d of the upper triangular matrix
        xy = np.ones((n, n))*np.linspace(0, n-1, n) 
        x = xy.T[self.triu_idx]
        y = xy[self.triu_idx]
        
        # extract the upper triangular matrix of data into a flat array      
        z = self.data[self.triu_idx]
        
        if self.remove_nan:
            # remove NaN's
            x = x[~self.nan_idx_triu]
            y = y[~self.nan_idx_triu]
            z = z[~self.nan_idx_triu]
        
        if self.scaler:
            scaler.fit(z.reshape(1, -1))
            z_scaled = scaler.transform(z.reshape(1, -1))
            z = z_scaled[0, :]
        #if log_scale == 2:
        #    z = np.log2(z)
        #elif log_scale == 10:
        #    z = np.log10(z)
            

        # create matrix with x, y, z as columns
        X = np.array([x, y, z]).T

        return X
    
    def scipy_dist(self, metric, col1, col2, **kwargs):
        """
        Computes the distance matrix using scipy.pdist with a given metric.
        
        Args:
            metric {'row_col_dist'}:
                Which metric to use when calculating the dissimilarity matrix
        Return:
            distmat (array):
                Distance matrix on square form
        """
        metric = "_" + metric
        if hasattr(self, metric):
            # call metric function
            metric = getattr(self, metric)
            
        # calculate the pairwise distances between data point and convert it
        # to square distance matrix
        distmat = squareform(pdist(self.X[:, col1:col2], metric=metric, **kwargs))

        return distmat
    
    def sklearn_dist(self, metric, scaler=None, **kwargs):
        """
        Computes the distance matrix using sklearn.pairwise_distances.
        
        Args:
            metric {'euclidean', 'test_euclidean'}:
                Which metric to use when calculating the dissimilarity matrix
        Returns:
            distmat (array):
                Distance matrix on square form
        """
        metric = "_" + metric
        if hasattr(self, metric):
            # call metric function
            metric = getattr(self, metric)

        # calculate the pairwise distances between data point and convert it
        # to square distance matrix
        distmat = pairwise_distances(self.X, metric=metric, force_all_finite=False, 
                                     n_jobs=-1, **kwargs)
        
        return distmat
        
    def _interactions_dist(self, u, v):
        """Calculates the difference in number of interactions.
        
        Calculates the difference in number of interactions given by
            dist = |u_z - v_z|
        where z is the number of interactions.
        
        Args:
            u (ndarray):
                Array containing one data point.
            v (ndarray):
                Array containing one data point.
        Result:
            dist (float):
                Difference in number of interactions between u and v.
        """
        dist = abs(u[2] - v[2])
        return dist
        
    def _diagonal_dist(self, u, v):
        """Calculates the scaled difference in distance from the main diagonal.
        
        Calculates the scaled  difference in distance from the main diagonal 
        between u and v given by
            dist = | |u_i - v_i| - |u_j - v_j| | / n,
        where i is the row index, j is the column index and n is the maximum 
        distance from the main diagonal. The maximum diagonal distance equals
        the length of the matrix when square.

        Args:
            u (ndarray):
                Array containing one data point.
            v (ndarray):
                Array containing one data point.
        Result:
            dist (float):
                Difference in distance from the main diagonal between u and v.
        """
        d1 = abs(u[0] - u[1])
        d2 = abs(v[0] - v[1])
        
        dist = abs(d1 - d2)/self.n
        return dist
    
    def _manhattan(self, u, v):
        """Calculates the scaled manahattan distance between u and v.
        
        Calculates the scaled manhattan distance between u and v given by
            dist = (|u_i - v_i| + |u_j - v_j|)/2n,
        where i is the row index,j is the column index and 2n is the maximum 
        manhattan distance.
        
        Args:
            u (ndarray):
                Array containing one data point.
            v (ndarray):
                Array containing one data point.
        Result:
            dist (float):
                Manahattan distance between u and v.
        """
        dist = (abs(u[0] - v[0]) + abs(u[1] - v[1]))/(2*self.n)
        return dist
    
    def _diag_3d_dist(self, u, v):
        d1 = np.sqrt((u[0] - u[1])**2/(2*self.n) + u[2]**2)
        d2 = np.sqrt((v[0] - v[1])**2/(2*self.n) + v[2]**2)
        
        dist = abs(d1 - d2)
          
    


if __name__ == '__main__':
    test = np.array([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])
    distmat = DissimilarityMatrix(test)
    print(distmat.scipy_dist('minkowski_', log_scale=None, k=2))
