import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import os
from pysam.libcvcf import defaultdict
import seaborn as sns


import plot_utils as pltut
import publish_plot as pubplt

def parse_args():
    parser = argparse.ArgumentParser(description='Plot vs. XXX Results')
    parser.add_argument('--continuous', type=str, required=True,
                        help='continuous output data file')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')
    return parser.parse_args()

def read_continuous(continuous):
    f = open(continuous, 'r')
    pval_data = []
    snps_data = []
    time_data = []
    header = f.readline()

    for line in f:
        L = line.strip().split(',')
        gene = L[0]
        pval = int(L[1])
        snps = int(L[2])
        time = float(L[3])

        pval_data.append(pval)
        snps_data.append(snps)
        time_data.append(time)

    return pval_data, snps_data, time_data

def plot_histogram(data,
                   label,
                   out):
    fig, ax = plt.subplots(1, 1, figsize=(10, 5), dpi=300)

    ax.hist(data, bins=100, color='cadetblue')

    # title
    # plt.suptitle('Histogram of {}'.format(label))
    ax.set_xlabel(label)
    ax.set_ylabel('Frequency')
    # ax.set_yscale('log')

    # add text: number of files
    ax.text(0.5, 0.5, 'Number of files: {}'.format(len(data)/19181),
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes)


    plt.legend(frameon=False)

    # format
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.savefig(out)



def plot_scatter(counts_data, time_data,
                 label, out):

    # scatter plot 1
    fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=300)
    # fig = plt.figure(figsize=(10, 5), dpi=300)
    # ax1 = fig.add_gridspec(1, 1, left=0.1, right=0.55, top=0.9, bottom=0.1).subplots()
    x = time_data
    y = counts_data

    # pval = ax.scatter(x, y1, alpha=0.1,
    #                    color='cadetblue',
    #                    label='pval hits')
    # hexagonal binning plot
    hb = ax.hexbin(x, y, gridsize=100, cmap='viridis',
                   mincnt=1, bins='log')

    # add colorbar to plot
    cb = plt.colorbar(hb, ax=ax)
    cb.set_label('counts in bin')



    # snp = ax.scatter(x, y2, alpha=0.2,
    #                   color='yellowgreen',
    #                   label='snps hits')

    # # plot histogram of xxx times on right rotated
    # # left=0.1, right=0.75, top=.9, bottom=0.1
    # # 0.75: This is the left position of the inset axes, meaning the inset starts at 75% of the width of the figure from the left.
    # # 0: This is the bottom position of the inset axes, meaning it starts at the bottom of the figure (0%).
    # # 0.25: This is the width of the inset axes, meaning it occupies 25% of the figure's width.
    # # 1: This is the height of the inset axes, meaning it occupies the full height of the figure (100%).
    # print('plotting pval histograms...')
    # ax_histy1 = ax1.inset_axes([.4, 0, 0.45, 1], sharey=ax1)
    # ax_histy1.hist(y1, bins=100,
    #               orientation='horizontal',
    #               color='cadetblue',
    #               alpha=0.5)
    # ax_histy1.set_xticks([])
    # plt.setp(ax_histy1.get_yticklabels(), visible=False)
    # ax_histy1.set_xscale('log')
    # ax_histy1.set_xticks([10, 1000, 100000])
    # ax_histy1.spines['top'].set_visible(False)
    # ax_histy1.spines['right'].set_visible(False)
    # print('plotting snps histograms...')
    # ax_histy2 = ax1.inset_axes([0.25, 0, 0.85, 1], sharey=ax1)
    # ax_histy2.hist(y2, bins=100,
    #               orientation='horizontal',
    #               color='yellowgreen',
    #               alpha=0.5)
    # ax_histy2.set_xticks([])
    # plt.setp(ax_histy2.get_yticklabels(), visible=False)
    # ax_histy2.set_xscale('log')
    # ax_histy2.set_xticks([10, 1000, 100000])
    # ax_histy2.spines['top'].set_visible(False)
    # ax_histy2.spines['right'].set_visible(False)

    # title
    ax.set_xlabel('XXX decompression time')
    ax.set_ylabel('XXX ' + label)

    # plt.legend(frameon=False)

    # format
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.savefig(out)



def main():

    args = parse_args()
    continuous_file = args.continuous
    out = args.out

    print('reading file...')
    pval_data, snps_data, time_data = read_continuous(continuous_file)

    print('plotting pvals...')
    plot_scatter(pval_data, time_data, 'pval hits', out + 'continuous-pval-hex.png')
    print('plotting snps...')
    plot_scatter(snps_data, time_data, 'SNPs hits', out + 'continuous-snps-hex.png')

    # print('plotting times...')
    # plot_histogram(time_data,
    #                'Decompression time by gene (s)',
    #                out + 'continuous-times-hist.png')
    #
    # print('plotting pvals...')
    # plot_histogram(pval_data,
    #                'Number of pval hits by gene',
    #                out + 'continuous-pval-hist.png')
    #
    # print('plotting snps...')
    # plot_histogram(snps_data,
    #                'Total number of significant SNPs by gene',
    #                out + 'continuous-snps-hist.png')


if __name__ == '__main__':
    main()