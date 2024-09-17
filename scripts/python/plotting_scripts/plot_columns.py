import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

def plot_column_scatter(column_sizes,
                        column_decompression_times,
                        colors,
                        out_png):

    # one scatter plot of column compressed sizes vs decompression times for each block size
    # color = codec cocktail
    fig, ax = plt.subplots(4, 1, figsize=(15, 7), dpi=300)
    for i, block_size in enumerate(column_sizes):
        for codec in column_sizes[block_size]:
            ax[i].scatter(column_sizes[block_size][codec],
                          column_decompression_times[block_size][codec],
                          color=colors[codec],
                          s=20,
                          marker='o',
                          edgecolor='black',
                          alpha=0.5)
            ax[i].set_title(block_size)
            ax[i].set_xlabel('Column compressed size (bytes)')
            ax[i].set_ylabel('Decompression time (seconds)')

    # FORMATTING
    # legend, no frame
    # add custom legend for codecs
    codec_legend_elements = [Patch(facecolor=colors[codec], edgecolor='black', label=codec)
                            for codec in colors]
    ax[0].legend(handles=codec_legend_elements, loc='upper left', frameon=False)

    # remove spines
    for i in range(4):
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)

    # SAVE
    plt.savefig(out_png)
