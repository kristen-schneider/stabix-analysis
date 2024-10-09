import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pysam.libcvcf import defaultdict

import plot_utils as pltu

def parse_args():
    parser = argparse.ArgumentParser(description='Plot file sizes')
    parser.add_argument('--decomp', type=str, required=True,
                        help='dir with decompression times')
    parser.add_argument('--colors', type=str, required=True,
                        help='path to file which lists for each codec')
    parser.add_argument('--out', type=str, required=True,
                        help='Output directory for plots')
    return parser.parse_args()


def plot_column_by_data_type_scatter(column_sizes,
                                     column_decompression_times,
                                     block_sizes,
                                     colors,
                                     out):
    data_types = ['int', 'float', 'string']
    fig, ax = plt.subplots(len(block_sizes), len(data_types),
                           figsize=(23, 20), dpi=300,
                           sharex=False, sharey=False)
    codecs = []


    for i, block_size in enumerate(block_sizes):
        for j, data_type in enumerate(data_types):
            for codec in column_sizes[block_size]:
                if codec not in codecs:
                    codecs.append(codec)
                try:
                    ax[i, j].scatter(column_sizes[block_size][codec][data_type],
                                 column_decompression_times[block_size][codec][data_type],
                                 color=colors[codec],
                                 s=20,
                                 alpha=0.5)
                    # log scale y-axis
                    ax[i, j].set_yscale('log')
                except KeyError:
                    print(f'No data for {block_size} {codec} {data_type}')


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
                            for codec in codecs]
    ax[0, 0].legend(handles=codec_legend_elements, loc='upper left', frameon=False)

    # remove spines
    for i in range(len(block_sizes)):
        for j in range(len(data_types)):
            ax[i, j].spines['top'].set_visible(False)
            ax[i, j].spines['right'].set_visible(False)

    # FORMATTING
    # plt.tight_layout()

    out_name = out + 'column_decompression.png'
    # SAVE
    plt.savefig(out_name)

def main():
    args = parse_args()
    colors = pltu.read_colors(args.colors)

    file_names = ['continuous-103220-both_sexes']
    codecs = ['bz2', 'deflate', 'xz', 'zlib', 'zstd', 'combo-fbb']
    block_sizes = ['1000', '2000', '5000', '10000', 'map']

    decomp_col_times = defaultdict(dict)
    comp_col_sizes = defaultdict(dict)

    for file_name in file_names:
        for block_size in block_sizes:
            for codec in codecs:
                full_file_name = args.decomp + '/' + file_name + '_' + block_size + '_' + codec + '_column_decompression.csv'
                decomp_col_times, comp_col_sizes = pltu.read_column_decompression(full_file_name,
                                                                                  block_size,
                                                                                  decomp_col_times,
                                                                                  comp_col_sizes)
        plot_column_by_data_type_scatter(comp_col_sizes,
                                         decomp_col_times,
                                         block_sizes,
                                         colors,
                                         args.out)



if __name__ == '__main__':
    main()
