import numpy as np

def constant(x):
    y = 2*x
    return y

def change_with_noise(x):
    print(np.sum(x == 0)/(np.shape(x)[0]*np.shape(x)[1]))
    print(np.nanmax(x))
    y = x.copy()
    for i in range(np.shape(x)[0]):
        i = int(i)
        for j in range(np.shape(x)[1]):
            j = int(j)
            y[i, j] = np.nanmax([0, 2*x[i, j] + np.random.normal(0, 0.01) ])
    print(np.nanmax(y), np.nanmin(y), np.mean(y), np.sum(y == 0)/(np.shape(y)[0]*np.shape(y)[1]))
    return y