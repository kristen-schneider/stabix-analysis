import argparse
from collections import defaultdict
import os
from pysam.libcvcf import defaultdict


import plot_utils as pltut
import publish_plot as pubplt

def parse_args():
    parser = argparse.ArgumentParser(description='Write table of info for files')
    parser.add_argument('--data', type=str, required=True,
                        help='dir with data')
    parser.add_argument('--bed', type=str, required=True,
                        help='bed file genes')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for output table')
    return parser.parse_args()

def main():

    args = parse_args()
    data = args.data
    bed = args.bed
    out = args.out

    out_names = ['combo-xbb']
    block_sizes = ['2000']

    files = pltut.read_pvals(data + '/pvals.txt').keys()
    genes = pltut.read_bed_file(bed)

    tabix_out_dir = data + '/tabix_output/'
    data_dir = data + 'data/'

    # get tabix data

    compressed_files = defaultdict()
    all_gene_times = defaultdict()
    all_gene_records = defaultdict()
    all_gene_pval_hits = defaultdict()
    all_genomic_indexes = defaultdict()

    num_files = 0

    for f in files:
        tabix_gene_times, tabix_gene_records, tabix_gene_pval_hits = (
            pltut.read_genes(tabix_out_dir + f + '_tabix_output.txt', False))
        all_gene_times[f] = {'tabix': tabix_gene_times[f]}
        all_gene_records[f] = {'tabix': tabix_gene_records[f]}
        all_gene_pval_hits[f] = {'tabix': tabix_gene_pval_hits[f]}
        all_genomic_indexes[f] = {}
        for block_size in block_sizes:
            for name in out_names:
                base_name = f + '_' + block_size + '_' + name
                file_dir = data_dir + base_name + '/'
                try:
                    tsv = data_dir + f + '.tsv'
                    tsv_size = os.path.getsize(tsv)
                    bgz = data_dir + f + '.tsv.bgz'
                    bgz_size = os.path.getsize(bgz)
                    tbi = data_dir + f + '.tsv.bgz.tbi'
                    tbi_size = os.path.getsize(tbi)
                    xxx = file_dir + base_name + '.grlz'
                    xxx_size = os.path.getsize(xxx)
                    gen = file_dir + 'genomic.idx'
                    gen_size = os.path.getsize(gen)
                    pval = file_dir + 'pval.idx'
                    pval_size = os.path.getsize(pval)

                    genomic_index = pltut.read_genomic_index(gen)
                    all_genomic_indexes[f][block_size] = {name: genomic_index}

                    xxx_query = file_dir + f + '.query'
                    xxx_gene_times, xxx_gene_records, xxx_gene_pval_hits = \
                        pltut.read_genes(xxx_query, True)
                    try:
                        all_gene_times[f][block_size] = {name: xxx_gene_times[f]}
                        all_gene_records[f][block_size] = {name: xxx_gene_records[f]}
                        all_gene_pval_hits[f][block_size] = {name: xxx_gene_pval_hits[f]}
                    except KeyError:
                        all_gene_times[f] = {block_size: {name: xxx_gene_times[f]}}
                        all_gene_records[f] = {block_size: {name: xxx_gene_records[f]}}
                        all_gene_pval_hits[f] = {block_size: {name: xxx_gene_pval_hits[f]}}

                    compressed_files[f] = {'tsv': tsv_size,
                                           'bgz': bgz_size, 'tbi': tbi_size,
                                           'xxx': xxx_size, 'gen': gen_size, 'pval': pval_size}
                except FileNotFoundError:
                    print('File not found: ' + base_name)

                num_files += 1

        # if num_files == 3:
        #     break

    pltut.write_all_files_table(all_gene_times,
                                all_gene_records,
                                all_gene_pval_hits,
                                all_genomic_indexes,
                                compressed_files,
                                out_names,
                                block_sizes,
                                genes,
                                out + 'table.csv')

if __name__ == '__main__':
    main()