import cooler

import numpy as np

class HiCInformation:
    def __init__(self, filename):
        """
        Args:
            filename (string):  
                Path to a multi-resolution coolers .mcool file
        """
        self.filename = filename
        
    def available_resolutions(self):
        """Prints available resolutions"""
        resolutions_str = cooler.fileops.list_coolers(self.filename)
        
        resolutions_int = np.zeros_like(resolutions_str, dtype=int)
        for idx, elem in enumerate(resolutions_str):
            resolutions_int[idx] = int(elem[13:])
            
        sorted_resolutions = np.sort(resolutions_int)
        print('Available resolutions:')
        print(*sorted_resolutions, sep="\n")