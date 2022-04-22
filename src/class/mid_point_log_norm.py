import numpy as np

from  matplotlib.colors import LogNorm

class MidPointLogNorm(LogNorm):
    def __init__(self, vmin=None, vmax=None, midpoint=None, log=np.log2, clip=False):
        LogNorm.__init__(self,vmin=vmin, vmax=vmax, clip=clip)
        self.midpoint=midpoint
        self.log = log
        
        
    def __call__(self, value, clip=None):
        result, is_scalar = self.process_value(value)
        self.autoscale_None(result)
        x, y = [self.log(self.vmin), self.log(self.midpoint), self.log(self.vmax)], [0, 0.5, 1]
        return np.ma.array(np.interp(self.log(result), x, y), mask=result.mask, copy=False)