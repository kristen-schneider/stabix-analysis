import argparse
import tabix_utils
import os
import sys


num_genes = 100
p_value_threshold = 5e-8
p_value_idx = 9

def parse_args():
    parser = argparse.ArgumentParser(description='Sandbox script')
    parser.add_argument('--bed', type=str, help='gene bed file')
    parser.add_argument('--gwas', type=str, help='GWAS file (.bgz) and tabix included')
    parser.add_argument('--pval_threshold', type=str, help='less than or equal to this value')
    parser.add_argument('--pval_cols', type=str, help='file with p-value columns')
    parser.add_argument('--tabix_out', type=str, help='output directory for tabix results')

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

    gwas_file_size = os.path.getsize(gwas_file)

    output_tabix_time_file = args.tabix_out + '-tabix_query_times.txt'
    # print('Output file: {}'.format(output_tabix_time_file))
    time_file = open(output_tabix_time_file, 'a')
    time_file.truncate(0)
    time_file.write('GWAS file: {}'.format(gwas_basename) + ', Bytes: {}\n'.format(gwas_file_size))
    time_file.write('Gene,Time(sec)\n')

    output_tabix_query_file = args.tabix_out + '-tabix_query_results.txt'
    # print('Output file: {}'.format(output_tabix_query_file))
    query_file = open(output_tabix_query_file, 'a')
    query_file.truncate(0)

    # if empty, print error and exit
    if not p_value_idx_dict:
        print('Error: p-value index file is empty')
        sys.exit(1)

    for p in p_value_idx_dict:
        # print gwas file and p-value index and label to query file
        query_file.write('GWAS file: {}'.format(gwas_basename) +
                         ',index: {}'.format(p) +
                         ',title: {}\n'.format(p_value_idx_dict[p]))

        query_idx = 0
        total_sig_snps = 0
        total_sig_genes = 0
        # write gwas file and size to query file
        for gene in genes:
            query_file.write('Gene: {}\n'.format(gene))
            for chrom in genes[gene]:
                for start, end in genes[gene][chrom]:
                    records = tabix_utils.tabix_query_with_threshold(time_file,
                                                                     gene,
                                                                     chrom,
                                                                     start,
                                                                     end,
                                                                     gwas_file,
                                                                     p-1,
                                                                     p_value_threshold)
                    if records:
                        split_records = [record.split() for record in records]
                        for record in split_records:
                            query_file.write('\t'.join(record) + '\n')
                        total_sig_snps += len(records)
                        total_sig_genes += 1
            query_idx += 1

        query_file.write('END' +
                         ',Genes: {}'.format(total_sig_genes) +
                         ',SNPs: {}\n'.format(total_sig_snps))



if __name__ == '__main__':
    main()
