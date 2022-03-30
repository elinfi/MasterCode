import numpy as np

from  matplotlib.colors import LogNorm

class MidPointLogNorm(LogNorm):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        LogNorm.__init__(self,vmin=vmin, vmax=vmax, clip=clip)
        self.midpoint=midpoint
        
        
    def __call__(self, value, clip=None):
        result, is_scalar = self.process_value(value)
        self.autoscale_None(result)
        x, y = [np.log2(self.vmin), np.log2(self.midpoint), np.log2(self.vmax)], [0, 0.5, 1]
        return np.ma.array(np.interp(np.log2(result), x, y), mask=result.mask, copy=False)