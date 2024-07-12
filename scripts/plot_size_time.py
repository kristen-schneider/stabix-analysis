import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import utils


def parse_args():
    parser = argparse.ArgumentParser(description='Plot column sizes')
    parser.add_argument('--file_sizes', type=str, required=True,
                        help='gwas file sizes')
    parser.add_argument('--timing_dir', type=str, required=True,
                        help='Directory containing timing data')
    parser.add_argument('--colors', type=str, required=True,
                        help='Colors for each codec')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for plots')
    return parser.parse_args()


def plot_timing_vs_size(timing_data,
                        file_sizes,
                        colors,
                        output_dir):
    '''
    Plot decompression time vs compressed data size
    :param timing_data:
    :param file_sizes:
    :param colors:
    '''

    markers = {2000: 'o', 5000: 's', 10000: 'D', 20000: 'P'}

    # x-axis: decompression time
    # y-axis: compressed data size
    # one subplot for each block size

    fig, axs = plt.subplots(len(markers), 1, figsize=(10, 12))
    for i, block_size in enumerate(timing_data):
        color_list = []
        decompression_times = []
        compressed_sizes = []
        for codec in timing_data[block_size]:
            decompression_times.append(timing_data[block_size][codec]['decompression'])
            compressed_sizes.append(file_sizes['GCST90179150_buildGRCh37.tsv.grlz_' + codec + '-' + str(block_size)])
            color_list.append(colors[codec])

        axs[i].scatter(decompression_times, compressed_sizes,
                    color=color_list,
                    s=100,
                    marker=markers[block_size],
                    edgecolor='black',
                    alpha=0.7,
                    label=block_size)

        axs[i].set_xlabel('Decompression Time (s)')
        axs[i].set_ylabel('Compressed Data Size (bytes)')
        axs[i].set_title('Decompression Time vs Compressed Data Size' + '\n(block size: ' + str(block_size) + ')')

        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)

        # legend for codecs
        codec_legend_elements = [Patch(facecolor=colors[codec], edgecolor='black', label=codec)
                                for codec in colors]

        if i == 0:
            axs[i].legend(handles=codec_legend_elements, title='Codec', loc='upper left', frameon=False)

    plt.tight_layout()
    plt.savefig(output_dir + 'decompression_time_vs_size.png')
    plt.close()


def main():
    args = parse_args()

    # Read the file sizes
    file_sizes = utils.read_gwas_file_sizes(args.file_sizes)

    codecs = ['bz2', 'deflate', 'fpfvb', 'xz', 'zlib', 'zstd']
    block_sizes = [2000, 5000]
    colors = utils.read_colors(args.colors)
    block_data = {block_size: defaultdict(dict) for block_size in block_sizes}

    # Read the timing data
    for block_size in block_sizes:
        timing_file = args.timing_dir + 'times_' + str(block_size) + '.txt'
        timing_data = utils.read_timing_file(timing_file,
                                             codecs)
        block_data[block_size] = timing_data


    # plot timing vs size
    plot_timing_vs_size(block_data,
                        file_sizes,
                        colors,
                        args.output_dir)


if __name__ == '__main__':
    main()
