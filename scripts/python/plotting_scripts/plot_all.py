import argparse
import matplotlib.pyplot as plt
import os
from collections import defaultdict

from scripts.python.plotting_scripts import plot_file_sizes as pfs
from scripts.python.plotting_scripts import plot_columns as pc

BLOCK_SIZES = [10000, 15000, 20000, 'map']

def parse_args():
    parser = argparse.ArgumentParser(description='Plot file sizes')
    parser.add_argument('--file_sizes', type=str, required=True,
                        help='gwas file sizes')
    parser.add_argument('--compressed_sizes', type=str, required=True,
                        help='file with compressed column sizes')
    parser.add_argument('--decompression_times', type=str, required=True,
                        help='file with decompression times')
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
    kzip_genomic_idx_sizes = defaultdict(dict)
    kzip_pvalue_idx_sizes = defaultdict(dict)

    with open(file_sizes_file, 'r') as f:
        for line in f:
            L = line.strip().split(',')
            if 'gwas' in L[0]:
                gwas_file = L[1]
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
                    kzip_sizes[gwas_file][block_size][codec_cocktail] = size
                except KeyError:
                    try:
                        kzip_sizes[gwas_file][block_size] = {}
                        kzip_sizes[gwas_file][block_size][codec_cocktail] = size
                    except KeyError:
                        kzip_sizes[gwas_file] = {}
                        kzip_sizes[gwas_file][block_size] = {}
                        kzip_sizes[gwas_file][block_size][codec_cocktail] = size

            elif L[0] == 'kzip_gen_idx':
                block_size = L[1]
                codec_cocktail = L[2]
                size = int(L[3])
                try:
                    kzip_genomic_idx_sizes[gwas_file][block_size][codec_cocktail] = size
                except KeyError:
                    try:
                        kzip_genomic_idx_sizes[gwas_file][block_size] = {}
                        kzip_genomic_idx_sizes[gwas_file][block_size][codec_cocktail] = size
                    except KeyError:
                        kzip_genomic_idx_sizes[gwas_file] = {}
                        kzip_genomic_idx_sizes[gwas_file][block_size] = {}
                        kzip_genomic_idx_sizes[gwas_file][block_size][codec_cocktail] = size

            elif L[0] == 'kzip_pval_idx':
                block_size = L[1]
                codec_cocktail = L[2]
                size = int(L[3])
                try:
                    kzip_pvalue_idx_sizes[gwas_file][block_size][codec_cocktail] = size
                except KeyError:
                    try:
                        kzip_pvalue_idx_sizes[gwas_file][block_size] = {}
                        kzip_pvalue_idx_sizes[gwas_file][block_size][codec_cocktail] = size
                    except KeyError:
                        kzip_pvalue_idx_sizes[gwas_file] = {}
                        kzip_pvalue_idx_sizes[gwas_file][block_size] = {}
                        kzip_pvalue_idx_sizes[gwas_file][block_size][codec_cocktail] = size
    return gzip_sizes, bgzip_sizes, kzip_sizes, kzip_genomic_idx_sizes, kzip_pvalue_idx_sizes

def read_column_sizes(column_sizes_file):
    column_sizes = defaultdict(dict)
    with open(column_sizes_file, 'r') as f:
        # read header
        f.readline()
        for line in f:
            L = line.strip().split(',')
            if 'gwas' in L[0]:
                continue
            else:
                block_size = L[0]
                column_codec = L[1]
                data_type = L[2]
                try:
                    column_sizes[block_size][column_codec][data_type] = [int(x) for x in L[3:]]
                except KeyError:
                    try:
                        column_sizes[block_size][column_codec] = {}
                        column_sizes[block_size][column_codec][data_type] = [int(x) for x in L[3:]]
                    except KeyError:
                        column_sizes[block_size] = {}
                        column_sizes[block_size][column_codec] = {}
                        column_sizes[block_size][column_codec][data_type] = [int(x) for x in L[3:]]
    return column_sizes

def read_column_decompression_times(decompression_times_file):
    decompression_times = defaultdict(dict)
    with open(decompression_times_file, 'r') as f:
        # read header
        f.readline()
        for line in f:
            L = line.strip().split(',')
            if 'gwas' in L[0]:
                continue
            else:
                block_size = L[0]
                column_codec = L[1]
                data_type = L[2]
                try:
                    decompression_times[block_size][column_codec][data_type] = [float(x) for x in L[3:]]
                except KeyError:
                    try:
                        decompression_times[block_size][column_codec] = {}
                        decompression_times[block_size][column_codec][data_type] = [float(x) for x in L[3:]]
                    except KeyError:
                        decompression_times[block_size] = {}
                        decompression_times[block_size][column_codec] = {}
                        decompression_times[block_size][column_codec][data_type] = [float(x) for x in L[3:]]
    return decompression_times

def main():
    args = parse_args()
    colors = read_colors(args.colors)

    # COMPRESSED FILE SIZES
    (gzip_sizes,
     bgzip_sizes,
     kzip_sizes,
     kzip_genomic_idx_sizes,
     kzip_pvalue_idx_sizes)  = read_file_sizes(args.file_sizes)

    # pfs.plot_file_sizes_scatter(gzip_sizes,
    #                             bgzip_sizes,
    #                             kzip_sizes,
    #                             kzip_genomic_idx_sizes,
    #                             kzip_pvalue_idx_sizes,
    #                             colors,
    #                             args.output_dir + '/file_sizes_scatter.png')
    pfs.plot_file_sizes_violin(gzip_sizes,
                            bgzip_sizes,
                            kzip_sizes,
                            kzip_genomic_idx_sizes,
                            kzip_pvalue_idx_sizes,
                            colors,
                            args.output_dir + '/file_sizes_violin.png')



    # # COLUMN DECOMPRESSION SPEED BY COMPRESSED SIZE
    # column_sizes = read_column_sizes(args.compressed_sizes)
    # column_decompression_times = read_column_decompression_times(args.decompression_times)
    # pc.plot_column_by_block_scatter(column_sizes,
    #                                 column_decompression_times,
    #                                 colors,
    #                                 args.output_dir + '/column_block_size_scatter.png')
    #
    # pc.plot_column_by_data_type_scatter(column_sizes,
    #                                     column_decompression_times,
    #                                     BLOCK_SIZES,
    #                                     colors,
    #                                     args.output_dir + '/column_data_type_scatter.png')





if __name__ == '__main__':
    main()