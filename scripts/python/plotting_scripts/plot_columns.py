import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

def plot_column_by_block_scatter(column_sizes,
                                 column_decompression_times,
                                 colors,
                                 out_png):

    data_type_markers = {'int': 'o',
                         'float': 's',
                            'string': '^'}

    # one scatter plot of column compressed sizes vs decompression times for each block size
    # color = codec cocktail

    block_sizes = list(column_sizes.keys())
    codecs = list(column_sizes[block_sizes[0]].keys())
    data_types = list(column_sizes[block_sizes[0]][codecs[0]].keys())

    fig, ax = plt.subplots(4, 1, figsize=(10, 15), dpi=300)

    for i, block_size in enumerate(block_sizes):
        for codec in codecs:
            for data_type in data_types:
                ax[i].scatter(column_sizes[block_size][codec][data_type],
                              column_decompression_times[block_size][codec][data_type],
                              color=colors[codec],
                              s=20,
                              marker=data_type_markers[data_type],
                              alpha=0.5)
                ax[i].set_title('Block Size: ' + block_size, fontsize=16, fontweight='bold')
                ax[i].set_xlabel('Column compressed size\n(bytes)')
                ax[i].set_ylabel('Decompression time\n(microseconds)')
    # for i, block_size in enumerate(column_sizes):
    #     for codec in column_sizes[block_size]:
    #         for data_type in column_sizes[block_size][codec]:
    #             print(block_size, codec, data_type)
    #             ax[i].scatter(column_sizes[block_size][codec][data_type],
    #                           column_decompression_times[block_size][codec][data_type],
    #                           color=colors[codec],
    #                           s=20,
    #                           marker=data_type_markers[data_type],
    #                           alpha=0.5)
    #             ax[i].set_title('Block Size: ' + block_size, fontsize=16, fontweight='bold')
    #             ax[i].set_xlabel('Column compressed size\n(bytes)')
    #             ax[i].set_ylabel('Decompression time\n(microseconds)')

    # FORMATTING
    # legend, no frame
    # add custom legend for codecs
    codec_legend_elements = [Patch(facecolor=colors[codec], edgecolor='black', label=codec)
                            for codec in colors]
    # add markers for data types
    data_type_legend_elements = [Line2D([0], [0], marker=data_type_markers[data_type],
                                        color='w', markerfacecolor='black', label=data_type)
                                for data_type in data_type_markers]
    ax[0].legend(handles=codec_legend_elements + data_type_legend_elements, loc='upper left', frameon=False)
    # ax[0].legend(handles=codec_legend_elements, loc='upper left', frameon=False)

    # remove spines
    for i in range(4):
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)

    # FORMATTING
    plt.tight_layout()

    # SAVE
    plt.savefig(out_png)

def plot_column_by_data_type_scatter(column_sizes,
                                     column_decompression_times,
                                     block_sizes,
                                     colors,
                                     out_png):
    # rows by block size
    # columns by data type
    data_types = 3
    fig, ax = plt.subplots(len(block_sizes), data_types,
                           figsize=(15, 15), dpi=300,
                           sharex=False, sharey=False)

    block_sizes = ['10000', '15000', '20000', 'map']
    codecs = list(column_sizes[block_sizes[0]].keys())
    data_types = list(column_sizes[block_sizes[0]][codecs[0]].keys())

    for i, block_size in enumerate(block_sizes):
        for j, data_type in enumerate(data_types):
            for codec in codecs:
                ax[i, j].scatter(column_sizes[block_size][codec][data_type],
                                 column_decompression_times[block_size][codec][data_type],
                                 color=colors[codec],
                                 s=20,
                                 alpha=0.5)



    # LABELING
    # title columns by data type
    for j, data_type in enumerate(data_types):
        ax[0, j].set_title('Data Type: ' + data_type, fontsize=16, fontweight='bold')
       # add x-axis to bottom row
        ax[len(block_sizes) - 1, j].set_xlabel('Column compressed size\n(bytes)')

    # title rows by block size
    for i, block_size in enumerate(block_sizes):
        # add text to left of left column
        ax[i, 0].text(-0.3, 0.5, 'Block Size: ' + block_size, fontsize=16, fontweight='bold',
                      rotation=90, verticalalignment='center', horizontalalignment='center',
                      transform=ax[i, 0].transAxes)
        # add y-axis to left column
        ax[i, 0].set_ylabel('Decompression time\n(microseconds)')

    # FORMATTING
    # legend, no frame
    # add custom legend for codecs
    codec_legend_elements = [Patch(facecolor=colors[codec], edgecolor='black', label=codec)
                            for codec in colors]
    ax[0, 0].legend(handles=codec_legend_elements, loc='upper left', frameon=False)

    # remove spines
    for i in range(len(block_sizes)):
        for j in range(len(data_types)):
            ax[i, j].spines['top'].set_visible(False)
            ax[i, j].spines['right'].set_visible(False)

    # FORMATTING
    # plt.tight_layout()

    # SAVE
    plt.savefig(out_png)
