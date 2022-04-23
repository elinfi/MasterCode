import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import pretty_plotting as pplot

def plot_rand_idx(rand_idxes, region, labels,
                  s=16, m=18, l=20,
                  figsize=(10, 7.5),
                  path='Replicates',
                  save=False):
    # use seaborn style
    sns.set_theme('paper')
    sns.set_style('ticks')
    
    # set font sizes
    pplot.font_size(s, m, l)
    
    f, axs = plt.subplots(figsize=figsize,
                          nrows=1,
                          ncols=1,
                          sharex=True, sharey=False)
    
    # adjust subplots
    plt.subplots_adjust(top=0.95,
                        bottom=0.17,
                        left=0.15,
                        right=0.99)
    
    ax = axs
    
    final_rand_idxes = rand_idxes[1:, :-1]
    mask = np.zeros_like(final_rand_idxes)
    mask[np.triu_indices_from(mask, 1)] = True
    
    sns.heatmap(final_rand_idxes, 
                vmin=0, vmax=1, 
                linewidths=.1, 
                annot=True, square=False, mask=mask,
                xticklabels=labels[:-1],
                yticklabels=labels[1:],
                ax=ax,
                cbar_kws={'pad': 0.02})
    ax.set_title('Rand index', y=0)
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    
    if save:
        plt.savefig(f'/home/elinfi/MasterCode/img/statistics/{region}_rand_idx.pdf')