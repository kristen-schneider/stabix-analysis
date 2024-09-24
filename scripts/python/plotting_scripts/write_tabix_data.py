import argparse
from collections import defaultdict
from os.path import basename

import plot_utils as plot_utils


def parse_args():
    parser = argparse.ArgumentParser(description='Plot Tabix Results')
    parser.add_argument('--times', type=str, required=True,
                        help='tabix search times')
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


def read_tabix_times(tabix_times_file):
    # gwas file: times
    tabix_times = defaultdict()
    f = open(tabix_times_file, 'r')
    for line in f:
        if len(line.strip()) == 0:
            continue
        elif "GWAS file: " in line:
            # read header1 and get gwas file name
            header_1 = line.split(',')
            gwas_file = header_1[0].strip().split(':')[1].strip()
            tabix_times[gwas_file] = 0.0
            # read header2
            header_2 = f.readline().split(',')
        else:
            L = line.strip().split(',')
            gene = L[0]
            time = float(L[1])
            tabix_times[gwas_file] += time

    return tabix_times


def write_tabix_time_data(tabix_times,
                          tabix_results,
                          out):

    out_file_name = out + 'tabix_search_times.csv'
    out_file = open(out_file_name, 'w')
    out_file.write('GWAS_file,Sig_Genes,Sig_SNPs,Total_Time(s)\n')
    for gwas_file, time in tabix_times.items():
        file_name = gwas_file.replace('.tsv.bgz', '')
        try:
            genes, snps = tabix_results[file_name][0]
        except KeyError:
            genes = -1
            snps = -1
        out_file.write(file_name + ','
                       + str(genes) + ','
                       + str(snps) + ','
                       + str(time) + '\n')
    out_file.close()


def read_tabix_results(tabix_results_file):
    tabix_results = defaultdict()
    f = open(tabix_results_file, 'r')
    gene_records = 0
    for line in f:
        L = line.strip().split(',')
        if 'GWAS file: ' in L[0]:
            gwas_file = line.strip().split(',')[0].split(':')[1].replace('.tsv.bgz', '').strip()
            tabix_results[gwas_file] = []
        elif 'Gene:' in L[0]:
            gene_records = 0
        elif 'END'in L[0]:
            genes = int(L[1].split(':')[1].strip())
            snps = int(L[2].split(':')[1].strip())
            try:
                tabix_results[gwas_file].append((genes, snps))
            except KeyError:
                tabix_results[gwas_file] = [(genes, snps)]
        else:
            # reading tabix results
            # print(gwas_file, line)
            gene_records += 1

    return tabix_results


def write_tabix_results(tabix_results,
                        out):
    out_file_name = out + 'tabix_search_results.csv'
    out_file = open(out_file_name, 'w')
    out_file.write('GWAS File,Genes,SNPs\n')
    for gwas_file, results in tabix_results.items():
        for result in results:
            out_file.write(gwas_file + ','
                           + str(result[0]) + ','
                           + str(result[1]) + '\n')
    out_file.close()


def main():
    args = parse_args()

    gene_sizes = get_gene_size(args.bed)

    tabix_times = read_tabix_times(args.times)
    tabix_results = read_tabix_results(args.results)

    write_tabix_time_data(tabix_times,
                          tabix_results,
                          args.out)

    # write_tabix_results(tabix_results,
    #                     args.out)


if __name__ == '__main__':
    main()