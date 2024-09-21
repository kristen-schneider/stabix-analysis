import argparse
import matplotlib.pyplot as plt
import os
from collections import defaultdict

import plot_utils as pu

def parse_args():
    parser = argparse.ArgumentParser(description='Plot Tabix Results')
    parser.add_argument('--times', type=str, required=True,
                        help='tabix search times')
    parser.add_argument('--bed', type=str, required=True,
                        help='bed file with gene sizes')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')

    return parser.parse_args()

def plot_tabix_times(tabix_times,
                     out):
    # for each gwas file, plot a histogram of tabix search times
    num_gwas = len(tabix_times)
    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)
    for i, (gwas_file, times) in enumerate(tabix_times.items()):
        curr_times = times.values()

        ax.hist(curr_times, bins=100, alpha=0.5, label=gwas_file)
        ax.set_xlabel('Tabix search time per gene (s)')
        ax.set_ylabel('Frequency')
        ax.set_title('Tabix Search Times')

        # format plot
        ax.legend(frameon=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)


        fig.savefig(os.path.join(out, 'tabix_search_times.png'))
        plt.close(fig)

def plot_tabix_times_gene_size(tabix_times,
                               gene_sizes,
                               out):
    # for each gwas file, plot a scatter plot of tabix search times (X) vs gene size (Y)
    num_gwas = len(tabix_times)
    times_data = []
    sizes_data = []
    for gwas_file, times in tabix_times.items():
        for gene, time in times.items():
            times_data.append(time)
            sizes_data.append(gene_sizes[gene])

    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)
    ax.scatter(sizes_data, times_data, alpha=0.5)
    ax.set_xlabel('Gene Size (bp)')
    ax.set_ylabel('Tabix search time per gene (s)')
    ax.set_title('Tabix Search Times vs Gene Size')

    # format plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.savefig(os.path.join(out, 'tabix_search_times_gene_size.png'))
    plt.close(fig)



def main():
    args = parse_args()

    tabix_times = pu.read_tabix_times(args.times)
    gene_sizes = pu.get_gene_size(args.bed)

    plot_tabix_times(tabix_times,
                     args.out)
    plot_tabix_times_gene_size(tabix_times,
                               gene_sizes,
                               args.out)





if __name__ == '__main__':
    main()

