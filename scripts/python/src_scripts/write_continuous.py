import argparse
import os, sys
from collections import defaultdict

sys.path.append('scripts/python/plotting_scripts/')
import plot_utils as pltut

def parse_args():
    parser = argparse.ArgumentParser(description='Plot Tabix vs. XXX Results')
    parser.add_argument('--data', type=str, required=True,
                        help='dir with data')
    parser.add_argument('--bed', type=str, required=True,
                        help='bed file genes')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')
    return parser.parse_args()


def main():
    args = parse_args()

    data = args.data
    bed = args.bed
    out = args.out

    genes = pltut.read_bed_file(bed)

    # gene: [(pval_hits, num_records, time)]
    all_gene_info = {}
    # all_gene_records = {}
    # all_gene_pval_hits = {}

    for f in os.listdir(data):
        if 'continuous' in f:
            query_file = data + f + '/' + f.replace('_2000_combo-xzb','') + '.query'
            basename = query_file.split('/')[-1].replace('.query','')
            if 'kristen' in f:
                x = 1
            try:
                gene_times, gene_records, gene_pval_hits = pltut.read_genes(query_file, True)
                for gene in genes:
                    try:
                        all_gene_info[gene].append((gene_pval_hits[basename][gene],
                                                    gene_records[basename][gene],
                                                    gene_times[basename][gene]))
                    except KeyError:
                        try:
                            all_gene_info[gene] = [(gene_pval_hits[basename][gene],
                                                    gene_records[basename][gene],
                                                    gene_times[basename][gene])]
                        except KeyError:
                            try:
                                all_gene_info[gene].append((0, 0, 0))
                            except KeyError:
                                all_gene_info[gene] = [(0, 0, 0)]
                print('done: ', f)
            except FileNotFoundError:
                continue
            except NotADirectoryError:
                continue

    # write gene info to file
    out_file = out + 'gene_info.csv'
    header = "Gene,Pval Hits,Num Records,Time(s)\n"
    with open(out_file, 'w') as f:
        f.write(header)
        for gene in all_gene_info:
            for info in all_gene_info[gene]:
                f.write("{},{},{},{}\n".format(gene, info[0], info[1], info[2]))





if __name__ == '__main__':
    main()