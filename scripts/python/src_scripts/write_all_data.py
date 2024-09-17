## write all data which can be used for plotting
## Variables:
## block sizes
## codecs
## data types

import argparse
import os
from collections import defaultdict

import numpy as np
import sys

# GLOBAL VARIABLES
CODEC_COCKTAILS = ['bz2', 'deflate', 'xz']
BLOCK_SIZES = [10000, 15000, 20000, 'map']

# argumnnet parser
def parse_args():
    parser = argparse.ArgumentParser(description='Writing GWAS Compression Data')
    parser.add_argument('--root', type=str, required=True,
                        help='Directory containing GWAS files and subdirectories.')
    parser.add_argument('--gwas', type=str, required=True,
                        help='Text file containing GWAS file basenames.')
    parser.add_argument('--comp', type=str, required=True,
                        help='Directory containing compressed timimgs.')
    parser.add_argument('--decomp', type=str, required=True,
                        help='Directory containing decompressed timings.')
    parser.add_argument('--out', type=str, required=True,
                        help='Output directory for written data.')
    return parser.parse_args()

# read gwas file names
def read_gwas_files(gwas_file):
    '''
    Read in gwas file names from a text file.
    '''
    gwas_files = []
    with open(gwas_file, 'r') as f:
        for line in f:
            gwas_files.append(line.strip())
    return gwas_files

def get_file_sizes(gwas_files,
                   root,
                   block_sizes,
                   codecs):
    # get file sizes for bgzip and gzip
    # gwas file: block_size: size
    gzip_sizes = {}
    bgzip_sizes = {}

    for gwas_file in gwas_files:
        # bgzip
        bgzip_file = os.path.join(root,
                                  gwas_file + '.tsv.bgz')
        bgzip_file_tbxi = os.path.join(root,
                                        gwas_file + '.tsv.bgz.tbi')
        bgzip_sizes[gwas_file] = os.path.getsize(bgzip_file) + os.path.getsize(bgzip_file_tbxi)
        # gzip
        gzip_file = os.path.join(root,
                                 gwas_file + '.tsv.gz')
        gzip_sizes[gwas_file] = os.path.getsize(gzip_file)

    # get file size for kzip
    # gwas file: block_size: codec: size
    kzip_sizes = defaultdict(dict)
    for gwas_file in gwas_files:
        for block_size in BLOCK_SIZES:
            for codec_cocktail in CODEC_COCKTAILS:
                kzip_file = os.path.join(root,
                                         gwas_file + '_' + str(block_size) + '_' + codec_cocktail,
                                         gwas_file + '_' + str(block_size) + '_' + codec_cocktail + '.grlz')
                kzip_file_gen_idx = os.path.join(root,
                                                 gwas_file + '_' + str(block_size) + '_' + codec_cocktail,
                                                 'genomic.idx')
                try:
                    kzip_sizes[gwas_file][block_size][codec_cocktail] = (
                            os.path.getsize(kzip_file) + os.path.getsize(kzip_file_gen_idx))
                except KeyError:
                    try:
                        kzip_sizes[gwas_file][block_size] = {}
                        kzip_sizes[gwas_file][block_size][codec_cocktail] = (
                            os.path.getsize(kzip_file)) + os.path.getsize(kzip_file_gen_idx)
                    except KeyError:
                        kzip_sizes[gwas_file] = {}
                        kzip_sizes[gwas_file][block_size] = {}
                        kzip_sizes[gwas_file][block_size][codec_cocktail] = (
                            os.path.getsize(kzip_file)) + os.path.getsize(kzip_file_gen_idx)

    return gzip_sizes, bgzip_sizes, kzip_sizes

def get_decompression_results(decomp_dir):
    '''
    Read in decompression results from a directory.
    '''
    # gwas file: block_size: codec: size OR time
    col_comp_sizes = {}
    col_decomp_times = {}

    for file in os.listdir(decomp_dir):
        gwas_basename = file.split('.')[0]
        gwas_basename = file.split('_')[0] + '-' + file.split('_')[1]

        if file.endswith('column_decompression.csv'):
            # get block size and codec
            block_size, codec_cocktail = file.split('_')[2], file.split('_')[3]
            # read in column compressed sizes
            col_sizes, col_times = read_decompression_results(os.path.join(decomp_dir, file))
            try:
                col_comp_sizes[gwas_basename][block_size][codec_cocktail] = col_sizes
            except KeyError:
                try:
                    col_comp_sizes[gwas_basename][block_size] = {}
                    col_comp_sizes[gwas_basename][block_size][codec_cocktail] = col_sizes
                except KeyError:
                    col_comp_sizes[gwas_basename] = {}
                    col_comp_sizes[gwas_basename][block_size] = {}
                    col_comp_sizes[gwas_basename][block_size][codec_cocktail] = col_sizes
            try:
                col_decomp_times[gwas_basename][block_size][codec_cocktail] = col_times
            except KeyError:
                try:
                    col_decomp_times[gwas_basename][block_size] = {}
                    col_decomp_times[gwas_basename][block_size][codec_cocktail] = col_times
                except KeyError:
                    col_decomp_times[gwas_basename] = {}
                    col_decomp_times[gwas_basename][block_size] = {}
                    col_decomp_times[gwas_basename][block_size][codec_cocktail] = col_times

    return col_comp_sizes, col_decomp_times

def read_decompression_results(decomp_file):
    # codec: [col_size, col_size, ...]
    col_sizes = {}
    # codec: [decomp_time, decomp_time, ...]
    col_decomp_times = {}

    # open file and read in data, split by comma
    # block_idx, col_idx, comp_time, col_size, codec
    with open(decomp_file, 'r') as f:
        # read header
        header = f.readline()
        data = f.readlines()
    # get the column sizes and decompression times
    col_sizes = {}
    col_decomp_times = {}
    for line in data:
        block_idx, col_idx, comp_time, col_size, codec = line.strip().split(',')
        # remove microseconds unit from compression time
        comp_time = float(comp_time.replace('μs', ''))
        try:
            col_sizes[codec].append(col_size)
        except KeyError:
            col_sizes[codec] = [col_size]
        try:
            col_decomp_times[codec].append(comp_time)
        except KeyError:
            col_decomp_times[codec] = [comp_time]
    return col_sizes, col_decomp_times


def main():

    args = parse_args()

    # read gwas file names
    gwas_files = []
    with open(args.gwas, 'r') as f:
        for line in f:
            gwas_files.append(line.strip())

    # get file sizes for bgzip, gzip, and kzip
    gzip_sizes, bgzip_sizes, kzip_sizes = get_file_sizes(gwas_files,
                                                         args.root,
                                                         BLOCK_SIZES,
                                                         CODEC_COCKTAILS)

    # write file sizes to output directory
    # gwas file: bgzip_size, gzip_size, kzip_size
    with open(os.path.join(args.out, 'file_sizes.csv'), 'w') as f:
        # gwas file
        # gzip size
        # bgzip size
        # block size: codec_cocktail: kzip size
        for gwas_file in gwas_files:
            f.write('gwas file,' + f'{gwas_file}' + '\n')
            f.write('gzip,' + str(gzip_sizes[gwas_file]) + '\n')
            f.write('bgzip,' + str(bgzip_sizes[gwas_file]) + '\n')
            for block_size in BLOCK_SIZES:
                for codec_cocktail in CODEC_COCKTAILS:
                    try:
                        f.write('kzip,' +
                                str(block_size) + ',' +
                                codec_cocktail + ',' +
                                str(kzip_sizes[gwas_file][block_size][codec_cocktail]) + '\n')
                    except KeyError:
                        f.write('NA,')
            f.write('\n')


    # get size and timing for kzip columns
    col_sizes, col_decomp_times = get_decompression_results(os.path.join(args.root, args.decomp))
    x = 'pause'

    # get overall compression times for kzip
    # get overall decompression times for kzip



if __name__ == '__main__':
    main()
