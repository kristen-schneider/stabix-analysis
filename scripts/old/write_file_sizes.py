import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data_dir', help='directory gwas files', required=True)
    parser.add_argument('-g', '--UKBB_gwas_files.txt', help='list of gwas files to plot', required=True)
    return parser.parse_args()

def read_gwas_files_names(gwas_names_file):
    '''
    Read the names of the gwas files from a file
    :param gwas_names_file: text file containing the names of the gwas files to include in plot
    :return: list of gwas file names
    '''
    gwas_files = []
    with open(gwas_names_file, 'r') as f:
        for line in f:
            gwas_files.append(line.strip())
    return gwas_files

def get_sizes_of_files(data_dir):
    '''
    Get the sizes of the files in the data directory
    :param data_dir: directory containing the gwas files (compressed and index files)
    :return: dictionary of file names and their sizes in bytes
    '''
    files = os.listdir(data_dir)
    file_sizes = {}
    for file in files:
        file_size_bytes = os.path.getsize(os.path.join(data_dir, file))
        # # convert bytes to MB
        # file_size_mb = file_size_bytes / 1024 / 1024
        # file_sizes[file] = file_size_mb
        file_sizes[file] = file_size_bytes
    return file_sizes

def write_file_sizes(gwas_file_names,
                     file_sizes,
                     output_file):
    '''
    Write the file sizes to a file
    :param gwas_file_names: list of gwas file names
    :param file_sizes: dictionary of file names and their sizes in bytes
    :param output_file: output file to write the file sizes
    :return: None
    '''
    with open(output_file, 'w') as f:
        f.write('File, Size (Bytes)\n')
        for file, size in file_sizes.items():
            for gwas_file in gwas_file_names:
                if gwas_file in file:
                    f.write(f'{file}, {size}\n')
    f.close()


def main():
    args = parse_args()
    data_dir = args.data_dir
    gwas_files = args.gwas_files

    block_sizes = [2000, 5000, 10000, 20000, 'map']
    codecs = ['bz2', 'deflate', 'fpfvb', 'zlib', 'zstd', 'xz']

    gwas_file_names = read_gwas_files_names(gwas_files)
    file_sizes = get_sizes_of_files(data_dir)

    output_file = os.path.join(data_dir, 'gwas_file_sizes.csv')
    write_file_sizes(gwas_file_names, file_sizes, output_file)



if __name__ == '__main__':
    main()