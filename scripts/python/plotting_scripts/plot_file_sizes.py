import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

def plot_file_sizes(gzip_sizes,
                    bgzip_sizes,
                    kzip_sizes,
                    codec_colors,
                    out_png):
    # get average gzip size
    avg_gzip_size = np.mean(gzip_sizes)
    # get average bgzip size
    avg_bgzip_size = np.mean(bgzip_sizes)

    # get average kzip sizes for each block size and codec
    avg_kzip_sizes = {}
    for block_size in kzip_sizes:
        avg_kzip_sizes[block_size] = {}
        for codec in kzip_sizes[block_size]:
            avg_kzip_sizes[block_size][codec] = np.mean(kzip_sizes[block_size][codec])

    # plot the average sizes
    fig, ax = plt.subplots(1, 1, figsize=(15, 7), dpi=300)
    # y = size in bytes
    # x = block size
    # color = codec cocktail

    # plot the average gzip size in red horizontal line
    plt.axhline(y=avg_gzip_size,
                color='black', linestyle='solid', linewidth=3,
                label='gzip')

    # plot the average bgzip size in blue horizontal line
    plt.axhline(y=avg_bgzip_size,
                color='black', linestyle='dashed', linewidth=3,
                label='bgzip')

    # single point for kzip codec cocktail sizes
    for block_size in avg_kzip_sizes:
        for codec in avg_kzip_sizes[block_size]:
            plt.scatter(block_size, avg_kzip_sizes[block_size][codec],
                        color=codec_colors[codec],
                        s=200,
                        marker='o',
                        edgecolor='black',
                        alpha=0.5)


    # FORMATTING
    # legend, no frame
    # add custom legend for codecs
    codec_legend_elements = [Patch(facecolor=codec_colors[codec], edgecolor='black', label=codec)
                            for codec in codec_colors]
    # add bgzip dashed line to legend
    codec_legend_elements.insert(0, Line2D([0], [0],
                                           color='black', linestyle='dashed', linewidth=3,
                                           label='bgzip'))
    # add gzip solid line to legend
    codec_legend_elements.insert(1, Line2D([0], [0],
                                           color='black', linestyle='solid', linewidth=3,
                                           label='gzip'))

    plt.legend(handles=codec_legend_elements,
               title='Compression Method', loc='best', frameon=False,
               fontsize=12, title_fontsize=14)

    # label the y-axis
    plt.ylabel('Compressed file size\n(bytes)', fontsize=14)

    # label the x-axis
    x_ticks = range(len(kzip_sizes.keys()))
    x_tick_labels = [f'{block_size}' for block_size in kzip_sizes.keys()]
    plt.xticks(x_ticks, x_tick_labels, fontsize=12)
    plt.xlabel('Block size\n(number of records in block)', fontsize=14)

    # remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # SAVE
    plt.savefig(out_png)
