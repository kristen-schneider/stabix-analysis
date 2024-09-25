import argparse

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
    parser.add_argument('--pval_threshold', type=str, help='less than or equal to this value')
    parser.add_argument('--pval_cols', type=str, help='file with p-value columns')
    parser.add_argument('--tabix_out', type=str, help='test_output directory for tabix results')

    return parser.parse_args()

def main():

    args = parse_args()

    bed_file = args.bed
    genes = tabix_utils.get_genes(bed_file)
    gene_count = len(genes)
    # print('Number of genes: {}'.format(gene_count))

    gwas_file = args.gwas
    gwas_basename = os.path.basename(gwas_file)
    p_value_idx_dict = tabix_utils.read_pval_index(args.pval_cols)
    p_value_threshold = float(args.pval_threshold)

    output_tabix_query_file = args.tabix_out + '-tabix_query_results.txt'

    out_file = open(output_tabix_query_file, 'a')
    out_file.truncate(0)

    # if empty, print error and exit
    if not p_value_idx_dict:
        print('Error: p-value index file is empty')
        sys.exit(1)

    for p in p_value_idx_dict:
        # header for results file
        out_file.write('GWAS file: {}'.format(gwas_basename.replace('.tsv.bgz','')) +
                         ',index: {}'.format(p-1) +
                         ',title: {}\n'.format(p_value_idx_dict[p]))

        total_sig_snps = 0
        total_sig_genes = 0
        # write gwas file and size to query file
        gene_i = 0
        for gene in genes:
            records = []
            if gene == "CSF2RA":
                x = 0
            gene_results = []
            for chrom in genes[gene]:
                for start, end in genes[gene][chrom]:
                    records, time = tabix_utils.tabix_query_with_threshold(chrom,
                                                                           start,
                                                                           end,
                                                                           gwas_file,
                                                                           p-1,
                                                                           p_value_threshold)
                    gene_results.append((records, time))
                    if records:
                        total_sig_snps += len(records)
                        total_sig_genes += 1
            # write gene and time to query file
            gene_total_time = sum([time for _, time in gene_results])
            out_file.write('Gene: {}'.format(gene) + ',time: {}\n'.format(gene_total_time))
            # write gene results to query file
            # print(gene_i, gene, len(records), total_sig_snps)
            for records, _ in gene_results:
                split_records = [r.split() for r in records]
                for record in split_records :
                    out_file.write('\t'.join(record) + '\n')

            gene_i += 1

        out_file.write('END' +
                         ',GWAS file: {}'.format(gwas_basename.replace('.tsv.bgz','')) +
                         ',index: {}'.format(p) +
                         ',title: {}'.format(p_value_idx_dict[p]) +
                         ',Genes: {}'.format(total_sig_genes) +
                         ',SNPs: {}\n'.format(total_sig_snps))



if __name__ == '__main__':
    main()
