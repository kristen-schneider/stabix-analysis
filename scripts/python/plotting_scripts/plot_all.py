import argparse
import matplotlib.pyplot as plt
import os
from collections import defaultdict

from scripts.python.plotting_scripts import plot_file_sizes as pfs


def parse_args():
    parser = argparse.ArgumentParser(description='Plot file sizes')
    parser.add_argument('--file_sizes', type=str, required=True,
                        help='gwas file sizes')
    parser.add_argument('--colors', type=str, required=True,
                        help='path to file which lists for each codec')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for plots')
    return parser.parse_args()

def read_colors(colors_file):
    # codec, color
    colors = {}
    with open(colors_file, 'r') as f:
        for line in f:
            L = line.strip().split(',')
            colors[L[0]] = L[1]
    return colors


def read_file_sizes(file_sizes_file):
    gzip_sizes = []
    bgzip_sizes = []
    kzip_sizes = defaultdict(dict)

    with open(file_sizes_file, 'r') as f:
        for line in f:
            L = line.strip().split(',')
            if 'gwas' in L[0]:
                continue
            elif L[0] == 'gzip':
                size = int(L[1])
                gzip_sizes.append(size)
            elif L[0] == 'bgzip':
                size = int(L[1])
                bgzip_sizes.append(size)
            elif L[0] == 'kzip':
                block_size = L[1]
                codec_cocktail = L[2]
                size = int(L[3])
                try:
                    kzip_sizes[block_size][codec_cocktail] = size
                except KeyError:
                    kzip_sizes[block_size] = {}
                    kzip_sizes[block_size][codec_cocktail] = size
    return gzip_sizes, bgzip_sizes, kzip_sizes




def main():
    args = parse_args()

    gzip_sizes, bgzip_sizes, kzip_sizes = read_file_sizes(args.file_sizes)

    colors = read_colors(args.colors)

    pfs.plot_file_sizes(gzip_sizes,
                        bgzip_sizes,
                        kzip_sizes,
                        colors,
                        args.output_dir + '/file_sizes.png')





if __name__ == '__main__':
    main()