import kmedoids

import numpy as np
import numpy.ma as ma

from numba import jit
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import MinMaxScaler



class DissimilarityMatrix:
    def __init__(self, data):
        """Creates a dissimilarity matrix.

        Creates a dissimilarity based on a Hi-C data matrix with a given
        distance metric.

        Args:
            data (ndarray):
                Square symmetric matrix containing the Hi-C data to cluster
            metric {'euclidean', 'test_euclidean'}:
                Which metric to use when calculating the dissimilarity matrix
        Attributes:
            n (int):
                Length of data matrix
        """
        self.data = data
        self.n = data.shape[0]
        
    def scalar_projection(self, a, b):
        """Calculates the scalar projection of a on b.
        
        Args:
            a (ndarray):
                Vector coordinates of a
            b (ndarray):
                Vector coordinates of unit vector b
        Return:
            scalar_proj (ndarrray):
                Scalar projection of a on b
        """
        # scalar_proj = (np.dot(a, b)/np.linalg.norm(b))
        scalar_proj = a - np.dot(a, b)*b
        
        return scalar_proj
    def row_col_dist(self, u, v, w1, w2):
        """Calculates the weighted distance between u and v (row - col).
        
        Calculates the weighted distance between u and v
            w1*|z_u - z_v| + w2*| |x_u - x_v| - |y_u - y_v| |,
        where z is the number of interactions, x is the row index and y is the  
        column index. The sum of the weights equals 1.

        Args:
            u (ndarray):
                Array containing one data point.
            v (ndarray):
                Array containing one data point.
            w1 (float):
                Interaction weight.
            w2 (float):
                Diagonal distance weight.
        Result:
            dist (float):
                Weighted distance between points u and v.
        """
        d1 = abs(u[0] - u[1])
        d2 = abs(v[0] - v[1])
        d3 = abs(u[2] - v[2])
        
        dist = w1*d3 + w2*abs(d1 - d2)/self.n
        return dist
        
    def projection_dist(self, u, v, diag_unit, w1, w2):
        """Calculates the weighted distance between u and v (projection).
        
        Calculates the weighted distance between u and v
            w1*|z_u - z_v| + w2*|proj_u - proj_v|,
        where z is the number of interactions and proj_u and proj_v is the 
        projection of u and v on the diagonal, respectively. The sum of the
        weights equals 1.

        Args:
            u (ndarray):
                Array containing one point.
            v (ndarray):
                Array containing one point.
            w1 (float):
                Interaction weight.
            w2 (float):
                Diagonal distance weight.
        Result:
            dist (float):
                Weighted distance between points u and v.
        """       
        # longest distance from the diagonal
        #max_dist_diag =  np.linalg.norm(diag)/2
        
        # size of diagonal clusters
        #size_diag_clusters = max_dist_diag/medoids
        
        max_diag = self.n/np.sqrt(2)
        
        # the maximum distance from the diagonal
        proj_u = np.linalg.norm(self.scalar_projection(u[:2], diag_unit))
        proj_v = np.linalg.norm(self.scalar_projection(v[:2], diag_unit))
        
        d3 = abs(u[2] - v[2])
        dist = w1*d3 + w2*abs(proj_u - proj_v)/max_diag
        
        return dist

    def minkowski_(self, u, v, k):
        """Calculates the minkowski distance of order k between point u and v.

        Args:
            u (ndarray):
                Array containing one point.
            v (ndarray):
                Array containing one point.
        Result:
            dist (float):
                Euclidean distance between the two points u and v
        """
        u_v = u - v
        dist = np.linalg.norm(u_v, ord=k)
        return dist

    def transform_data_nan(self, scaler=None):
        """Extract the upper triangular of the data matrix.

        Return:
            X (ndarray):
                n x 3 array containing the bin1 and bin2 in the first two
                columns and the number of interactions in the last column
        """
        #assert log_scale in [None, 'log2', 'log10', 'MinMax'], 'Only None, 2 or 10 are valid log scales'
            
        n = self.n
        self.triu_idx = np.triu_indices(n)
        self.nan_idx = np.isnan(self.data)
        self.triu_nan_idx = self.nan_idx[self.triu_idx]

        # create 2d coordinate system in 1d of the upper triangular matrix
        xy = np.ones((n, n))*np.linspace(0, n-1, n) 
        x = xy.T[self.triu_idx]
        x = x[~self.triu_nan_idx]
        y = xy[self.triu_idx]
        y = y[~self.triu_nan_idx]
        
        # extract the upper triangular matrix of data into a flat array      
        z = self.data[self.triu_idx]
        z = z[~self.triu_nan_idx]
        if scaler:
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
 
    def transform_data(self, scaler=None):
        """Extract the upper triangular of the data matrix.

        Return:
            X (ndarray):
                n x 3 array containing the bin1 and bin2 in the first two
                columns and the number of interactions in the last column
        """
        #assert log_scale in [None, 'log2', 'log10', 'MinMax'], 'Only None, 2 or 10 are valid log scales'
            
        n = self.n
        triu_idx = np.triu_indices(n)

        # create 2d coordinate system in 1d of the upper triangular matrix
        xy = np.ones((n, n))*np.linspace(0, n-1, n)
        x = xy.T[triu_idx]
        y = xy[triu_idx]
        
        # extract the upper triangular matrix of data into a flat array      
        z = self.data[triu_idx]
        if scaler:
            scaler.fit(z.reshape(1, -1))
            z_scaled = scaler.transform(z.reshape(1, -1))
            z = z_scaled[0, :]
        #if log_scale == 2:
        #    z = np.log2(z)
        #elif log_scale == 10:
        #    z = np.log10(z)
            

        # create matrix with x, y, z as columns
        X = np.array([x, y, z]).T
        
        masked_X = ma.masked_invalid(X)

        return masked_X
    
    def tranform_data2(self, scaler=None):
        n = self.n
        triu_idx = np.triu_indices(n)
        masked_data = ma.masked_invalid(self.data)
        
        masked_triu_idx = masked_data(triu_idx).mask
        
        # create 2d coordinate system in 1d of the upper triangular matrix
        xy = np.ones((n, n))*np.linspace(0, n-1, n)
        x = xy.T[triu_idx]
        y = xy[triu_idx]
        
        # extract the upper triangular matrix of data into a flat array      
        z = self.data[triu_idx]
        if scaler:
            scaler.fit(z.reshape(1, -1))
            z_scaled = scaler.transform(z.reshape(1, -1))
            z = z_scaled[0, :]
     
        # create matrix with x, y, z as columns
        X = ma.array([x, y, z]).T
        X.mask = np.stack((masked_triu_idx, masked_triu_idx, masked_triu_idx), axis=1)
        
        

    def scipy_dist(self, metric, scaler=None, remove_nan=True, **kwargs):
        """
        Computes the distance matrix using scipy.pdist.
        
        Args:
            metric {'euclidean', 'test_euclidean'}:
                Which metric to use when calculating the dissimilarity matrix
        Return:
            distmat (array):
                Distance matrix on square form
        """
        # transform data
        if remove_nan:
            X = self.transform_data_nan(scaler)
        else:
            X = self.transform_data(scaler)
        
        if hasattr(self, metric):
            # call metric function
            metric = getattr(self, metric)
            
        # calculate the pairwise distances between data point and convert it
        # to square distance matrix
        distmat = squareform(pdist(X, metric=metric, **kwargs))

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
        # transform data
        X = self.transform_data(scaler)        
        
        if hasattr(self, metric):
            # call metric function
            metric = getattr(self, metric)

        # calculate the pairwise distances between data point and convert it
        # to square distance matrix
        distmat = pairwise_distances(X, metric=metric, force_all_finite=False, n_jobs=-1, **kwargs)
        
        return distmat


if __name__ == '__main__':
    test = np.array([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])
    distmat = DissimilarityMatrix(test)
    print(distmat.scipy_dist('minkowski_', log_scale=None, k=2))
