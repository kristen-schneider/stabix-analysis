import argparse

import matplotlib.pyplot as plt
import os
from collections import defaultdict

import plot_utils as plot_utils

def parse_args():
    parser = argparse.ArgumentParser(description='Plot Tabix Results')
    parser.add_argument('--times', type=str, required=True,
                        help='tabix search times')
    parser.add_argument('--results', type=str, required=True,
                        help='tabix search results')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')
    return parser.parse_args()

def plot_tabix_times_hist(tabix_times,
                          num_genes,
                          pval,
                          out):

    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)
    all_file_times = []
    for gwas_file, (genes, snps, total_time) in tabix_times.items():
        all_file_times.append(total_time)

    ax.hist(all_file_times, bins=10, alpha=0.5, density=True)

    ax.set_xlabel('Tabix Search Time (s)')
    ax.set_ylabel('Frequency')
    ax.set_title('Tabix Search Times', fontweight='bold')

    # format plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # add text box about the data
    num_gwas_files = len(tabix_times)
    text = 'Number of GWAS Files: {}\n'.format(num_gwas_files)
    text += 'Number of Genes: {}\n'.format(num_genes)
    text += 'P-value Threshold: {}\n'.format(pval)
    text += 'Mean Time: {:.2f}s\n'.format(sum(all_file_times) / len(all_file_times))
    ax.text(0.95, 0.95, text, verticalalignment='top', horizontalalignment='right',
            transform=ax.transAxes, fontsize=8,
            bbox=dict(facecolor='none', alpha=0.5, edgecolor='none'))

    plt.tight_layout()
    fig.savefig(out + 'tabix_search_times_hist.png')


def plot_tabix_times_by_num_hits(tabix_times,
                                 num_genes,
                                 pval,
                                 out):

    # plot violoin plot of number of Gene Hits (X) vs Tabix Search times (Y) for all files

    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)
    data = {}
    for gwas_file, (genes, snps, total_time) in tabix_times.items():
        try:
            data[genes].append(total_time)
        except KeyError:
            data[genes] = [total_time]

    ax.violinplot(data.values(), showmeans=True, showmedians=True)
    ax.set_xlabel('Number of Gene Hits\n(p-value < ' + str(pval) + ')')
    ax.set_ylabel('Tabix Search Time (s)')
    ax.set_title('Tabix Search Times by\nNumber of Gene Hits', fontweight='bold')

    # format plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # x-axis ticks
    ax.set_xticks(range(1, 4))
    ax.set_xticklabels(data.keys())

    # add text box about the data
    num_gwas_files = len(tabix_times)
    text = 'Number of GWAS Files: {}\n'.format(num_gwas_files)
    text += 'Number of Genes: {}\n'.format(num_genes)
    text += 'P-value Threshold: {}\n'.format(pval)
    ax.text(0.80, 0.45, text, verticalalignment='top', horizontalalignment='right',
            transform=ax.transAxes, fontsize=8, bbox=dict(facecolor='none', alpha=0.5, edgecolor='none'))

    plt.tight_layout()
    fig.savefig(out + 'tabix_search_times_by_gene_hits.png')


def read_tabix_results_data(tabix_results_file):
    header = None
    tabix_results = defaultdict()
    f = open(tabix_results_file, 'r')
    for line in f:
        if header == None:
            header = line.strip().split(',')
            continue
        L = line.strip().split(',')
        gwas_file = L[0]
        gene = int(L[1])
        snps = int(L[2])
        tabix_results[gwas_file] = (gene, snps)

    return tabix_results

def plot_tabix_results_hist(tabix_results,
                            out):
    # plot one histogram of the number of genes and one histogram of the number of snps for each file
    fig, ax = plt.subplots(2, 1, figsize=(10, 4), dpi=300)

    gene_data = []
    snp_data = []

    for gwas_file, (gene, snps) in tabix_results.items():
        gene_data.append(gene)
        snp_data.append(snps)

    ax[0].hist(gene_data, bins=5, alpha=0.5, density=True)
    ax[0].set_xlabel('Number of Significant Genes')
    ax[0].set_ylabel('Frequency')
    ax[0].set_xlim(0, 5)

    ax[1].hist(snp_data, bins=5, alpha=0.5, density=True)
    ax[1].set_xlabel('Number of Significant SNPs')
    ax[1].set_ylabel('Frequency')
    ax[1].set_xlim(0, 5)

    fig.suptitle('Tabix Results')

    # format plot
    for a in ax:
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)

    fig.savefig(out + 'tabix_results_hist.png')

def main():
    args = parse_args()
    num_genes = 20386
    pval = 5e-08

    print('Reading tabix times data...')
    tabix_times = plot_utils.read_tabix_times_data(args.times)
    print('...plotting times histograms...')
    plot_tabix_times_hist(tabix_times,
                          num_genes,
                          pval,
                          args.out)
    print('...plotting times scatter plot...')
    plot_tabix_times_by_num_hits(tabix_times,
                                 num_genes,
                                 pval,
                                 args.out)
    # print('...plotting times scatter plot...')
    # plot_tabix_times_by_gene_size(tabix_times,
    #                               args.out)
    #
    # print('Reading tabix results data...')
    # tabix_results = read_tabix_results_data(args.results)
    # print('...plotting results histograms...')
    # plot_tabix_results_hist(tabix_results,
    #                         args.out)


if __name__ == '__main__':
    main()

