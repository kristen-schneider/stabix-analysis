import argparse
import tabix_utils


num_genes = 100
p_value_threshold = 5e-8
p_value_idx = 9

def parse_args():
    parser = argparse.ArgumentParser(description='Sandbox script')
    parser.add_argument('--bed', type=str, help='gene bed file')
    parser.add_argument('--gwas', type=str, help='GWAS file (.bgz) and tabix included')
    parser.add_argument('--pval_index', type=int, help='index of the p-value column')
    parser.add_argument('--pval_threshold', type=str, help='less than or equal to this value')
    parser.add_argument('--tabix_out', type=str, help='output directory for tabix results')

    return parser.parse_args()

def main():

    args = parse_args()

    bed_file = args.bed
    genes = tabix_utils.get_genes(bed_file)
    gene_count = len(genes)
    print('Number of genes: {}'.format(gene_count))


    gwas_file = args.gwas
    p_value_idx = int(args.pval_index)
    p_value_threshold = float(args.pval_threshold)

    output_tabix_time_file = args.tabix_out + 'tabix_query_time.txt'
    time_file = open(output_tabix_time_file, 'a')
    time_file.truncate(0)
    time_file.write('Gene,Time(sec)\n')

    output_tabix_query_file = args.tabix_out + 'tabix_query_results.txt'
    query_file = open(output_tabix_query_file, 'a')
    query_file.truncate(0)

    query_idx = 0
    for gene in genes:
        print('Gene: {}'.format(gene) + ' ' + str(query_idx))
        query_file.write('Gene: {}\n'.format(gene))
        for chrom in genes[gene]:
            for start, end in genes[gene][chrom]:
                records = tabix_utils.tabix_query_with_threshold(time_file,
                                                                 gene,
                                                                 chrom,
                                                                 start,
                                                                 end,
                                                                 gwas_file,
                                                                 p_value_idx,
                                                                 p_value_threshold)
                if records:
                    for record in records:
                        query_file.write('\t'.join(record) + '\n')

        query_idx += 1


if __name__ == '__main__':
    main()
