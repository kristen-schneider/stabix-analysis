import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


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

def get_col_size_data(col_sizes_dict):
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

def plot_col_sizes(col_size_data,
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

    # create new output file for each block size
    for block_size in block_sizes:
        out_file = os.path.join(output_dir, f'col_sizes_{block_size}.png')

        # create a subplot for each data type
        fig, ax = plt.subplots(4, 1,
                               figsize=(10, 10), dpi=300)
                               # sharex=True, sharey=True)
        for i, data_type in enumerate(data_types):
            ax[i].set_title(data_type)
            block_data = []
            colors = []
            for codec in codecs:
                if codec not in col_size_data[data_type][block_size]:
                    continue
                colors.append(codec_colors[codec])
                block_data.append(col_size_data[data_type][block_size][codec])

            sns.kdeplot(data=block_data, ax=ax[i],
                        fill=True, common_norm=True, palette=colors)

            # legend for codecs
            legend_elements = [Line2D([0], [0], color=colors[i], lw=4, label=codec)
                                 for i, codec in enumerate(codecs) if codec in col_size_data[data_type][block_size]]
            ax[i].legend(handles=legend_elements, loc='upper right')



        plt.tight_layout()
        plt.savefig(out_file)


def read_colors(colors_file):
    '''
    read comma delimited colors file and store in a dict
    :param colors_file:
    :return:
    '''
    colors = {}
    with open(colors_file, 'r') as f:
        for line in f:
            codec, color = line.strip().split(',')
            colors[codec] = color

    return colors

def main():
    args = parse_args()
    col_size_files = os.listdir(args.col_size_dir)
    col_size_files = [f for f in col_size_files if f.endswith('.csv')]

    col_sizes_dict = defaultdict(dict)

    for col_size_file in col_size_files:
        col_size_file_path = os.path.join(args.col_size_dir, col_size_file)
        col_sizes_dict = read_col_size_file(col_size_file_path,
                                       col_sizes_dict)

    col_size_data = get_col_size_data(col_sizes_dict)
    codec_colors = read_colors(args.colors)
    plot_col_sizes(col_size_data,
                   codec_colors,
                   args.output_dir)



if __name__ == '__main__':
    main()

