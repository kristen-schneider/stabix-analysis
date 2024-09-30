import argparse
import time


from numpy.distutils.system_info import xft_info

import tabix_utils
import os
import sys


# num_genes = 100
# p_value_threshold = 5e-8
# p_value_idx = 9

def parse_args():
    parser = argparse.ArgumentParser(description='Sandbox script')
    parser.add_argument('--bed', type=str, help='gene bed file')
    parser.add_argument('--gwas', type=str, help='GWAS file (.bgz) and tabix included')
    parser.add_argument('--pval_threshold', type=str, help='greater than or equal to this value')
    parser.add_argument('--pval_index', type=str, help='file with p-value columns for each gwas file')
    parser.add_argument('--out', type=str, help='test_output directory for tabix results')

    return parser.parse_args()

def main():

    args = parse_args()

    bed_file = args.bed
    genes = tabix_utils.get_genes(bed_file)

    gwas_file = args.gwas
    gwas_file_basename = os.path.basename(gwas_file).replace('.tsv.bgz', '')

    p_value_idx_dict = tabix_utils.read_pval_index(args.pval_index)
    p_value = p_value_idx_dict[gwas_file_basename]
    # if empty, print error and exit
    if not p_value_idx_dict:
        print('Error: no p-value indexes found in {} for gwas file: {}'.format(args.pval_cols, gwas_file))
        sys.exit(1)

    p_value_threshold = float(args.pval_threshold)

    output_tabix_query_file = os.path.join(args.out)
    out_file = open(output_tabix_query_file, 'a')
    out_file.truncate(0)
    out_file.write('GWAS file: {}\n'.format(gwas_file_basename))

    total_sig_snps = 0
    total_sig_genes = 0
    # write gwas file and size to query file
    for gene in genes:
        for chrom in genes[gene]:
            for start, end in genes[gene][chrom]:
                start_time = time.time()
                records = tabix_utils.tabix_query_with_threshold(chrom,
                                                                       start,
                                                                       end,
                                                                       gwas_file,
                                                                       p_value,
                                                                       p_value_threshold)
                if records:
                    total_sig_snps += len(records)
                    total_sig_genes += 1
                # write gene and time to query file
                out_file.write('Gene: {}\n'.format(gene))
                # write snps to query file
                for r in records:
                    out_file.write(r + '\n')
                end_time = time.time()
                tabix_time = end_time - start_time
                out_file.write('Gene: {},time: {}\n'.format(gene, tabix_time))


if __name__ == '__main__':
    main()
