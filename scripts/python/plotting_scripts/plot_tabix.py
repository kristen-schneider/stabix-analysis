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

def read_tabix_times_data(tabix_times_file):
    # Bin,Gene,Gene_Size,Tabix_Times
    tabix_data = defaultdict()
    header = None
    f = open(tabix_times_file, 'r')
    for line in f:
        if header == None:
            header = line.strip().split(',')
            continue
        L = line.strip().split(',')
        bin = L[0]
        gene = L[1]
        gene_size = int(L[2])
        tabix_times = [float(x) for x in L[3:]]

        try:
            tabix_data[bin][gene_size].extend(tabix_times)
        except KeyError:
            try:
                tabix_data[bin][gene_size] = tabix_times
            except KeyError:
                tabix_data[bin] = {gene_size: tabix_times}

    return tabix_data


def plot_tabix_times_hist(tabix_times,
                          out):
    bins = [0, 0.5e9, 1e9, 1.5e9, 2e9, 2.5e9, 3e9]
    bin_names = ['0-0.5GB', '0.5-1GB', '1-1.5GB', '1.5-2GB', '2-2.5GB', '2.5-3GB']


    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)

    for bin, gene_sizes in tabix_times.items():
        bin_times = []
        for gene_size, times in gene_sizes.items():
            bin_times.extend(times)
        ax.hist(bin_times, bins=50, alpha=0.5, label=bin, density=True)

    ax.set_xlabel('Tabix search time (s)')
    ax.set_ylabel('Density')
    ax.set_title('Tabix Search Times by Gene Size')

    # format plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(frameon=False)

    fig.savefig(out + 'tabix_search_times_hist.png')

def plot_tabix_times_by_gene_size(tabix_times,
                                  out):

    # plot a scatter plot of tabix search times (X) vs gene size (Y) for all files

    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)
    for bin, gene_sizes in tabix_times.items():
        x_data = []
        y_data = []
        for gene_size, times in gene_sizes.items():
            avg_time = sum(times) / len(times)
            x_data.append(avg_time)
            y_data.append(gene_size)

        ax.scatter(x_data, y_data, alpha=0.5, label=bin)


    ax.set_xlabel('Average Tabix Search Time')
    ax.set_ylabel('Gene Size')
    ax.set_title('Tabix Search Times by Gene Size')

    # format plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(title='File Size', frameon=False)

    fig.savefig(out + 'tabix_search_times_by_gene_size_scatter.png')

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

    print('Reading tabix times data...')
    tabix_times = read_tabix_times_data(args.times)
    print('...plotting times histograms...')
    plot_tabix_times_hist(tabix_times,
                          args.out)
    print('...plotting times scatter plot...')
    plot_tabix_times_by_gene_size(tabix_times,
                                  args.out)

    print('Reading tabix results data...')
    tabix_results = read_tabix_results_data(args.results)
    print('...plotting results histograms...')
    plot_tabix_results_hist(tabix_results,
                            args.out)


if __name__ == '__main__':
    main()

