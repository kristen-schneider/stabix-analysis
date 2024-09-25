import argparse
from collections import defaultdict
from os.path import basename

import plot_utils as plot_utils


def parse_args():
    parser = argparse.ArgumentParser(description='Plot Tabix Results')
    parser.add_argument('--results', type=str, required=True,
                        help='tabix search results')
    parser.add_argument('--bed', type=str, required=True,
                        help='bed file with gene sizes')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')

    return parser.parse_args()


def get_gene_size(bed_file):
    gene_size = dict()
    f = open(bed_file, 'r')
    for line in f:
        L = line.strip().split()
        gene = L[3]
        size = int(L[2]) - int(L[1])
        gene_size[gene] = size
    return gene_size

def read_tabix_results(tabix_results_file):
    tabix_results = defaultdict()
    f = open(tabix_results_file, 'r')
    gene_records = 0
    for line in f:
        L = line.strip().split(',')
        if 'GWAS file: ' in L[0]:
            gwas_file = line.strip().split(',')[0].split(':')[1].strip()
            gwas_total_time = 0
        elif 'Gene:' in L[0]:
            gene_records = 0
            gene_name, gene_time = line.split(',')
            gene_name = gene_name.split(':')[1].strip()
            gene_time = float(gene_time.split(':')[1].strip())
            gwas_total_time += gene_time

        elif 'END'in L[0]:
            gwas_file = line.strip().split(',')[1].split(':')[1].strip()
            p_index = int(L[2].split(':')[1].strip())
            p_title = L[3].split(':')[1].strip()
            genes = int(L[4].split(':')[1].strip())
            snps = int(L[5].split(':')[1].strip())
            try:
                tabix_results[gwas_file].append((genes, snps, gwas_total_time))
            except KeyError:
                tabix_results[gwas_file] = [(genes, snps, gwas_total_time)]
        else:
            # reading tabix results
            # print(gwas_file, line)
            gene_records += 1

    return tabix_results


def write_tabix_results(tabix_results,
                        out):
    out_file_name = out + 'tabix_search_results.csv'
    out_file = open(out_file_name, 'w')
    out_file.write('GWAS_File,Sig_Genes,Sig_SNPs,Total_Time\n')
    for gwas_file, results in tabix_results.items():
        for result in results:
            out_file.write(gwas_file + ','
                           + str(result[0]) + ','
                           + str(result[1]) + ','
                            + str(result[2]) + '\n')
    out_file.close()


def main():
    args = parse_args()

    gene_sizes = get_gene_size(args.bed)

    tabix_results = read_tabix_results(args.results)

    write_tabix_results(tabix_results,
                          args.out)

    # write_tabix_results(tabix_results,
    #                     args.out)


if __name__ == '__main__':
    main()