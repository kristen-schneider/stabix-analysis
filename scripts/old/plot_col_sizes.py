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
    parser.add_argument('--col_size_dir', type=str, required=True,
                        help='Directory containing column sizes')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for plots')
    parser.add_argument('--colors', type=str, required=True,
                        help='Colors for each codec')
    return parser.parse_args()

def read_col_size_file(file_path, col_sizes_dict):
    '''

    :param file_path:
    :param col_sizes_dict:
    :return:
    '''
    # dict structure = codec, block_size, col_idx: [col_sizes]
    codec = file_path.split('/')[-1].split('_')[0]
    block_size = file_path.split('/')[-1].split('_')[1]

    f = open(file_path, 'r')
    header = f.readline()
    header2 = f.readline()
    print(codec, block_size)

    for line in f:
        L = line.strip().split(',')
        block_idx = int(L[0])
        col_idx = int(L[1])
        col_size = int(L[2])

        if codec not in col_sizes_dict:
            col_sizes_dict[codec] = {}
        if block_size not in col_sizes_dict[codec]:
            col_sizes_dict[codec][block_size] = {}
        if col_idx not in col_sizes_dict[codec][block_size]:
            col_sizes_dict[codec][block_size][col_idx] = []

        if col_idx == 1 and col_size > 100:
            print(f'block_idx: {block_idx}, col_size: {col_size}')

        col_sizes_dict[codec][block_size][col_idx].append(col_size)

    return col_sizes_dict

def get_col_size_by_data_type(col_sizes_dict):
    '''
    '''

    data_types_dict = {0: 'string',
                       1: 'chromosome',
                       2: 'basepair',
                       3: 'string', 4: 'string',
                       5: 'float', 6: 'float', 7: 'float', 8: 'float', 9: 'float'}

    # format to be data_type: block_size: codec: [col_sizes]
    col_size_data = {'string': defaultdict(dict),
                     'chromosome': defaultdict(dict),
                     'basepair': defaultdict(dict),
                     'float': defaultdict(dict)}

    for codec in col_sizes_dict:
        for block_size in col_sizes_dict[codec]:
            for col_idx in col_sizes_dict[codec][block_size]:
                col_sizes = col_sizes_dict[codec][block_size][col_idx]
                data_type = data_types_dict[col_idx]
                if block_size not in col_size_data[data_type]:
                    col_size_data[data_type][block_size] = {}
                if codec not in col_size_data[data_type][block_size]:
                    col_size_data[data_type][block_size][codec] = []
                col_size_data[data_type][block_size][codec].extend(col_sizes)


    return col_size_data

def plot_col_size_by_data_type(col_size_data,
                               codec_colors,
                               output_dir):
    '''
    make a separate plot for each datatype
    :param col_size_data:
    :param output_dir:
    :return:
    '''

    data_types = ['string', 'chromosome', 'basepair', 'float']
    codecs = ['bz2', 'deflate', 'xz', 'zlib', 'zstd']
    block_sizes = ['2000', '5000', '10000', '20000']

    # create new test_output file for each block size
    for block_size in block_sizes:
        out_file = os.path.join(output_dir, f'col_sizes_{block_size}_nofpf.png')

        # create a subplot for each data type
        fig, ax = plt.subplots(4, 1,
                               figsize=(10, 10), dpi=300)
                               # sharex=True)

        for i, data_type in enumerate(data_types):
            title = f'{data_type.capitalize()} data'
            ax[i].set_title(title, fontsize=16, fontweight='bold')
            ax[i].set_ylabel('Density', fontsize=14)
            block_data = []
            colors = []
            for codec in codecs:
                if codec not in col_size_data[data_type][block_size]:
                    continue
                colors.append(codec_colors[codec])
                block_data.append(col_size_data[data_type][block_size][codec])

            sns.kdeplot(data=block_data, ax=ax[i],
                        fill=False, common_norm=True, palette=colors)

            # legend for codecs
            if i == 0:
                legend_elements = [Line2D([0], [0], color=colors[i], lw=4, label=codec)
                                     for i, codec in enumerate(codecs) if codec in col_size_data[data_type][block_size]]
                ax[i].legend(handles=legend_elements, loc='upper right', frameon=False)
            else:
                ax[i].legend().set_visible(False)

            ax[i].set_xlabel('Column Size (bytes)', fontsize=14)


            # remove spines
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)


        plt.tight_layout()
        plt.savefig(out_file)

def get_col_size_by_column(col_sizes_dict):
    '''

    '''
    col_size_data = {col_idx: defaultdict(dict) for col_idx in range(10)}

    for codec in col_sizes_dict:
        for block_size in col_sizes_dict[codec]:
            for col_idx in col_sizes_dict[codec][block_size]:
                col_sizes = col_sizes_dict[codec][block_size][col_idx]
                if codec not in col_size_data[col_idx]:
                    col_size_data[col_idx][codec] = []
                col_size_data[col_idx][codec].extend(col_sizes)

    return col_size_data

def plot_col_size_by_column(col_size_data,
                            codec_colors,
                            output_dir):

    codecs = ['bz2', 'deflate', 'fpfvb', 'xz', 'zlib', 'zstd']
    block_sizes = ['2000', '5000', '10000', '20000']

    # five throuh 9: effect_allele_frequency	beta	standard_error	z	p_value
    col_name_dict = {0: 'ID', 1: 'Chromosome', 2: 'Position',
                     3: 'Effect Allele', 4: 'Other Allele',
                     5: 'Effect Allele Frequency', 6: 'Beta',
                     7: 'Standard Error', 8: 'Z', 9: 'P Value'}

    # create new test_output file for each block size
    for block_size in block_sizes:
        out_file = os.path.join(output_dir, f'col_sizes_{block_size}_fpf_column.png')

        # create a subplot for each data type
        fig, ax = plt.subplots(5, 2,
                               figsize=(10, 12), dpi=300)
        # sharex=True)

        # 5 rows, 2 columns
        for i, col_idx in enumerate(range(10)):
            title = f'Column {col_idx}: {col_name_dict[col_idx]}'
            ax[i // 2, i % 2].set_title(title, fontsize=16, fontweight='bold')
            ax[i // 2, i % 2].set_ylabel('Density', fontsize=14)
            block_data = []
            colors = []
            for codec in codecs:
                if codec not in col_size_data[col_idx]:
                    continue
                colors.append(codec_colors[codec])
                block_data.append(col_size_data[col_idx][codec])

            sns.kdeplot(data=block_data, ax=ax[i // 2, i % 2],
                        fill=False, common_norm=True, palette=colors)

            # legend for codecs
            if i == 0:
                legend_elements = [Line2D([0], [0], color=colors[i], lw=4, label=codec)
                                   for i, codec in enumerate(codecs) if codec in col_size_data[col_idx]]
                ax[i // 2, i % 2].legend(handles=legend_elements, loc='upper right', frameon=False)
            else:
                ax[i // 2, i % 2].legend().set_visible(False)

            ax[i // 2, i % 2].set_xlabel('Column Size (bytes)', fontsize=14)

            # remove spines
            ax[i // 2, i % 2].spines['top'].set_visible(False)
            ax[i // 2, i % 2].spines['right'].set_visible(False)


        plt.tight_layout()
        plt.savefig(out_file)


def main():
    args = parse_args()
    col_size_files = os.listdir(args.col_size_dir)
    col_size_files = [f for f in col_size_files if f.endswith('.csv')]

    codec_colors = utils.read_colors(args.colors)

    col_sizes_dict = defaultdict(dict)
    for col_size_file in col_size_files:
        col_size_file_path = os.path.join(args.col_size_dir, col_size_file)
        col_sizes_dict = read_col_size_file(col_size_file_path,
                                       col_sizes_dict)

    # col_sizes_by_data_type = get_col_size_by_data_type(col_sizes_dict)
    # plot_col_size_by_data_type(col_sizes_by_data_type,
    #                            codec_colors,
    #                            args.output_dir)

    col_sizes_by_column = get_col_size_by_column(col_sizes_dict)
    plot_col_size_by_column(col_sizes_by_column,
                            codec_colors,
                            args.output_dir)



if __name__ == '__main__':
    main()

