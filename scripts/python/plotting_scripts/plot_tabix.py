import argparse

import matplotlib.pyplot as plt
import os
from collections import defaultdict

import plot_utils as plot_utils

def parse_args():
    parser = argparse.ArgumentParser(description='Plot Results')
    parser.add_argument('--results', type=str, required=True,
                        help='search results')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')
    parser.add_argument('--name', type=str, required=False,
                        help='tabix/sqlite/hdf5')
    return parser.parse_args()

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

def plot_tabix_times_hist(tabix_times,
                          num_genes,
                          pval,
                          out,
                          name='tabix'):

    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)
    all_file_times = []
    for gwas_file, (genes, snps, total_time) in tabix_times.items():
        all_file_times.append(total_time)

    ax.hist(all_file_times, bins=10, alpha=0.5, density=True)

    ax.set_xlabel(f'{name.capitalize()} Search Time (s)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{name.capitalize()} Search Times', fontweight='bold')

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
    fig.savefig(out + f'{name}_search_times_hist.png')


def plot_tabix_times_by_num_hits(tabix_times,
                                 num_genes,
                                 pval,
                                 out,
                                 name='tabix'):

    # plot scatter plot of number of Genes (X) vs Tabix Search times (Y) for all files
    # plot scatter plot of number of SNPs (X) vs Tabix Search times (Y) for all files


    fig, ax = plt.subplots(2, 1, figsize=(8, 4), dpi=300)

    for gwas_file, (genes, snps, total_time) in tabix_times.items():
        ax[0].scatter(genes, total_time, alpha=0.5)
        ax[1].scatter(snps, total_time, alpha=0.5)

    # format plots
    for a in ax:
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)


    ax[0].set_xlabel('Number of Significant Genes')
    ax[0].set_ylabel(f'{name.capitalize()} Search Time (s)')
    # ax[0].set_title('Tabix Search Times by\nNumber of Significant Genes', fontweight='bold')

    ax[1].set_xlabel('Number of Significant SNPs')
    ax[1].set_ylabel(f'{name.capitalize()}Search Time (s)')
    # ax[1].set_title('Tabix Search Times by\nNumber of Significant SNPs', fontweight='bold')

    # add text box about the data
    num_gwas_files = len(tabix_times)
    text = 'Number of GWAS Files: {}\n'.format(num_gwas_files)
    text += 'Number of Genes: {}\n'.format(num_genes)
    text += 'P-value Threshold: {}\n'.format(pval)
    # ax.text(0.80, 0.45, text, verticalalignment='top', horizontalalignment='right',
    #         transform=ax.transAxes, fontsize=8, bbox=dict(facecolor='none', alpha=0.5, edgecolor='none'))

    plt.tight_layout()
    fig.savefig(out + f'{name}_search_times_by_gene_hits.png')

def plot_tabix_hits_hist(tabix_results,
                         num_genes,
                         pval,
                         out):
    # plot one histogram of the number of genes and one histogram of the number of snps for each file
    fig, ax = plt.subplots(2, 1, figsize=(10, 4), dpi=300)

    gene_data = []
    snp_data = []

    for gwas_file, (gene, snps, time) in tabix_results.items():
        gene_data.append(gene)
        snp_data.append(snps)

    ax[0].hist(gene_data, bins=20, alpha=0.5, density=True)
    ax[0].set_xlabel('Number of Significant Genes')
    ax[0].set_ylabel('Frequency')
    # ax[0].set_xlim(0, 5)

    ax[1].hist(snp_data, bins=20, alpha=0.5, density=True)
    ax[1].set_xlabel('Number of Significant SNPs')
    ax[1].set_ylabel('Frequency')
    # ax[1].set_xlim(0, 5)

    fig.suptitle(f'{name.capitalize()} Hits')

    # format plot
    for a in ax:
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)

    # add text box about the data
    num_gwas_files = len(tabix_results)
    text = 'Number of GWAS Files: {}\n'.format(num_gwas_files)
    text += 'Number of Genes: {}\n'.format(num_genes)
    text += 'P-value Threshold: {}\n'.format(pval)
    # ax[0].text(0.95, 0.95, text, verticalalignment='top', horizontalalignment='right',

    plt.tight_layout()
    fig.savefig(out + f'{name}_results_hist.png')

def main():
    args = parse_args()
    num_genes = 20386
    pval = 7.3
    name = args.name or 'tabix'

    print(f'Reading {name} times data...')
    tabix_times = plot_utils.read_tabix_times_data(args.results)
    print('...plotting times histograms...')
    plot_tabix_times_hist(tabix_times,
                          num_genes,
                          pval,
                          args.out,
                          name)
    print('...plotting hits histogram...')
    plot_tabix_hits_hist(tabix_times,
                         num_genes,
                         pval,
                         args.out,
                         name)
    print('...plotting times scatter plot...')
    plot_tabix_times_by_num_hits(tabix_times,
                                 num_genes,
                                 pval,
                                 args.out,
                                 name)


if __name__ == '__main__':
    main()

