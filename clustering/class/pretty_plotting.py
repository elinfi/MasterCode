import simulation_functions as sim

from matplotlib.ticker import EngFormatter

def format_ticks(ax, x=True, y=True, rotate=True):
    bp_formatter = EngFormatter('b')
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)
    return None
      
def region2extent(region):
    chrom, start, end = sim.split_region(region)
    
    extent = (start, end, end, start)
    
    return extent
