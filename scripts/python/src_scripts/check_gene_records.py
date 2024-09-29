import argparse
import os

import utils as utils

def parse_args():
    parser = argparse.ArgumentParser(description='Plot New Results')
    parser.add_argument('--data', type=str, required=True,
                        help='times dir')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')
    return parser.parse_args()

def check_gene_records(tabix_gene_records,
                       lzr_gene_records):

    # assert that there are the same number of genes in both records
    assert len(tabix_gene_records) == len(lzr_gene_records)

    for gene in tabix_gene_records.keys():
        if lzr_gene_records[gene] != tabix_gene_records[gene]:
            print('Gene:', gene)
            print('--LZR:', lzr_gene_records[gene])
            print('--Tabix:', tabix_gene_records[gene])


def main():
    args = parse_args()

    data_dir = args.data
    out_dir = args.out

    file_names = ['continuous-103220-both_sexes']
    # file_names = ['continuous-103220-both_sexes', 'continuous-30130-both_sexes-irnt']
    codecs = ['bz2']
    block_sizes = [5000]

    for file in file_names:
        file_dir = os.path.join(data_dir, file + '_output')
        tabix_file = os.path.join(file_dir, 'tabix_output.txt')
        tabix_gene_times, tabix_gene_records = utils.read_genes(tabix_file)
        for codec in codecs:
            for block_size in block_sizes:
                file_name = codec + '_' + str(block_size) + '.txt'
                file_path = os.path.join(file_dir, file_name)

                if not os.path.exists(file_path):
                    print('File does not exist:', file_path)
                    continue

                print('Reading file:', file_path)
                lzr_gene_times, lzr_gene_records = utils.read_genes(file_path)

                for gwas_file in tabix_gene_times.keys():
                    check_gene_records(tabix_gene_records[gwas_file],
                                       lzr_gene_records[gwas_file])



if __name__ == '__main__':
    main()