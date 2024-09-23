import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import utils

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--file_names', help='gwas file names', required=True)
    parser.add_argument('-s', '--file_sizes', help='gwas file sizes', required=True)
    parser.add_argument('-c', '--colors', help='colors for codecs', required=True)
    parser.add_argument('-p', '--plot_dir', help='test_output directory to save plots', required=True)
    return parser.parse_args()



def plot_file_sizes(gwas_file_names,
                    file_sizes,
                    block_sizes,
                    codec_colors,
                    plot_dir):
    '''
    Genearte a bar plot of the file sizes for each gwas file
    :param gwas_file_names: list of gwas file names
    :param file_sizes: dictionary of file names and their sizes in bytes
    :param block_sizes: list of block sizes
    :param codec_colors: dict of codecs and their colors
    :param plot_dir: test_output directory to save the plots
    :return: none
    '''

    # for each file, save a new plot
    # bar plot of file sizes, x-axis = block sizes, y-axis = file sizes
    # for each block size, there are many codecs to plot

    x_tick_i = 0

    codecs = sorted(list(codec_colors.keys()))

    for gwas_file in gwas_file_names:
        output_file = os.path.join(plot_dir, f'{gwas_file}_file_sizes.png')
        fig, ax = plt.subplots(1, 1, figsize=(15, 7), dpi=300)
        # for each block size, plot the file sizes
        for block_size in block_sizes:
            for codec in codecs:
                try:
                    file_name = f'{gwas_file}.tsv.grlz_{codec}-{block_size}'
                    file_size = file_sizes[file_name]
                    plt.bar(x_tick_i, file_size, color=codec_colors[codec])
                    # add text label for file size with 2 decimal places
                    float_size = file_size / 10e7
                    plt.text(x_tick_i, file_size,
                             f'{float_size:.2f}',
                             ha='center', va='bottom')
                    x_tick_i += 1
                except KeyError:
                    print(f'{file_name} not found in file sizes')
            x_tick_i += 1


        # x tick labels are the block sizes used and only need to be plotted once for each block size
        x_ticks = range(0, len(block_sizes) * len(codecs) + len(block_sizes), len(codecs) + 1)
        x_ticks = [x + len(block_sizes)/2 for x in x_ticks]
        # remove x-tick line
        ax.tick_params(axis='x', which='both', bottom=False, top=False)
        plt.xticks(x_ticks, block_sizes, fontsize=10)

        # plot gzip file size as a horizontal line
        gzip_file_name = f'{gwas_file}.tsv.gz'
        gzip_file_size = file_sizes[gzip_file_name]
        plt.axhline(y=gzip_file_size, color='black', linestyle='--', label='gzip')

        # log scale for y-axis
        plt.yscale('log')

        # remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.xlabel('Block Size', fontsize=15)
        plt.ylabel('File Size (Bytes)', fontsize=15)
        plt.title(f'File Sizes for {gwas_file}', fontsize=25)

        # legend is the codecs and their colors
        legend = [plt.Rectangle((0, 0), 1, 1, fc=codec_colors[codec]) for codec in codecs]
        legend.append(plt.Line2D([0], [0], color='black', linestyle='--', label='gzip'))
        legend_label = codecs.copy()
        legend_label.append('gzip')
        plt.legend(legend, legend_label, loc='upper right')
        # plt.tight_layout()

        plt.savefig(output_file)
        plt.close()


def main():
    args = parse_args()
    file_names = args.file_names
    file_sizes = args.file_sizes
    colors = args.colors
    plot_dir = args.plot_dir

    block_sizes = [2000, 5000, 10000, 20000, 'map']
    codec_colors = utils.read_colors(colors)

    gwas_file_names = utils.read_gwas_file_names(file_names)
    gwas_file_sizes = utils.read_gwas_file_sizes(file_sizes)

    plot_file_sizes(gwas_file_names,
                    gwas_file_sizes,
                    block_sizes,
                    codec_colors,
                    plot_dir)



if __name__ == '__main__':
    main()