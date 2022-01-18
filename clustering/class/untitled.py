import numpy as np

class ContactMatrix:
    def __init__(self):
        pass
    
    def transpose_coords(self, row, col):
        if row > col:
            return (row - col, row)
        else:
            return (col - row, col)