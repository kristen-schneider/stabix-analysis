import argparse

import matplotlib.pyplot as plt
import os
from collections import defaultdict
import plot_utils as plot_utils

def parse_args():
    parser = argparse.ArgumentParser(description='Plot Tabix vs. XXX Results')
    parser.add_argument('--tabix_files', type=str, required=True,
                        help='tabix search times by file')
    parser.add_argument('--new_files', type=str, required=True,
                        help='new search times by file')
    parser.add_argument('--tabix_genes', type=str, required=True,
                        help='tabix search times by gene')
    parser.add_argument('--new_genes', type=str, required=True,
                        help='new search times by gene')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')
    parser.add_argument('--compressed_files_dir', type=str, required=True,
                        help='directory with compressed files')
    return parser.parse_args()

def plot_compare_times_genes(tabix_genes,
                             new_genes,
                             tabix_gene_hits,
                             out):

    fig, ax = plt.subplots(1, 1, figsize=(10, 10), dpi=300)
    tabix_times_x= []
    new_times_y = []
    gene_records = []

    for gwas_file, genes in tabix_genes.items():
        if gwas_file not in new_genes:
            continue
        new_genes = new_genes[gwas_file]
        for gene, time in genes.items():
            tabix_times_x.append(time)
            new_times_y.append(new_genes[gene])
            gene_records.append(tabix_gene_hits[gwas_file][gene])

    plot = ax.scatter(tabix_times_x, new_times_y,
               c=gene_records, cmap='Paired',
               alpha=0.7)
    ax.set_ylim(0, 0.2)
    ax.set_xlim(0, 0.2)

    ax.set_xlabel('Tabix Search Time (s)')
    ax.set_ylabel('New Search Time (s)')
    ax.set_title('Tabix vs. New Search Times By Gene', fontweight='bold')

    # format
    cbar = plt.colorbar(plot)
    cbar.set_label('Number of Hits')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    # save
    plt.savefig(os.path.join(out, 'tabix_vs_new_search_times_genes.png'))



def plot_compare_times_files(tabix_times,
                             new_times,
                             num_genes,
                             pval,
                             out):
    # tabix times are in seconds (python time.time())
    # new times are in microseconds (c++ chrono::high_resolution_clock::now())

    # scatter: x = tabix time, y = new time, color = number of hits
    fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=300)

    for gwas_file, (genes, snps, tabix_time) in tabix_times.items():
        if gwas_file not in new_times:
            continue
        new_time = new_times[gwas_file]
        num_hits = genes
        ax.scatter(tabix_time, new_time, alpha=0.7)

    ax.set_xlabel('Tabix Search Time (s)')
    ax.set_ylabel('New Search Time (s)')
    ax.set_title('Tabix vs. New Search Times', fontweight='bold')

    # add text box about the data
    num_gwas_files = len(tabix_times)
    text = 'Number of GWAS Files: {}\n'.format(num_gwas_files)
    text += 'Number of Genes: {}\n'.format(num_genes)
    text += 'P-value Threshold: {}\n'.format(pval)
    ax.text(0.55, 0.85, text, verticalalignment='top', horizontalalignment='right',
            transform=ax.transAxes, fontsize=8,
            bbox=dict(facecolor='none', alpha=0.5, edgecolor='none'))

    # # legend
    # handles = [plt.Line2D([0], [0], marker='o', color='w',
    #                       markerfacecolor=hits_colors[i], label=i, markersize=8)
    #             for i in range(6)]
    #
    # ax.legend(handles=handles, title='Number of Hits', frameon=False, fontsize=8)

    # format
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    # save
    plt.savefig(os.path.join(out, 'tabix_vs_new_search_times.png'))


def plot_file_sizes(file_sizes,
                    out):

    # for each file, plot bar graph for each file size
    # x axis = file type
    # y axis = file size (GB)
    # color = file type
    # bar graph

    color_by_file_type = {'tsv': 'black',
                            'bgz': 'green',
                            'tbi': 'orange',
                            'new': 'blue',
                            'gen': 'purple',
                            'pval': 'red'}

    fig, ax = plt.subplots(1, 1, figsize=(10, 5), dpi=300)
    ax.set_title('File Sizes', fontsize=16)
    ax.set_ylabel('File Size (GB)', fontsize=14)

    x_points = range(len(file_sizes.keys()))
    width = 0.3

    full_files = ['tsv', 'bgz', 'new']
    for i, file_name in enumerate(file_sizes.keys()):
        x_offset = 0
        for j, file_type in enumerate(full_files):
            x_offset += width
            if file_type == '.bgz':
                ax.bar(i + x_offset, file_sizes[file_name][file_type],
                       color=color_by_file_type[file_type],
                       label=file_type, alpha=0.7, width=width)
                # stack '.tbi' on top
                ax.bar(i + x_offset, file_sizes[file_name]['tbi'],
                       color=color_by_file_type['tbi'],
                       label='tbi', alpha=0.7, width=width, bottom=file_sizes[file_name][file_type])
            elif file_type == 'new':
                ax.bar(i + x_offset, file_sizes[file_name][file_type],
                       color=color_by_file_type[file_type],
                       label=file_type, alpha=0.7, width=width)
                # stack 'gen' and 'pval' on top
                ax.bar(i + x_offset, file_sizes[file_name]['gen'],
                       color=color_by_file_type['gen'],
                       label='gen', alpha=0.7, width=width, bottom=file_sizes[file_name][file_type])
                ax.bar(i + x_offset, file_sizes[file_name]['pval'],
                       color=color_by_file_type['pval'],
                       label='pval', alpha=0.7, width=width, bottom=file_sizes[file_name][file_type] + file_sizes[file_name]['gen'])
            else:
                ax.bar(i + x_offset, file_sizes[file_name][file_type],
                       color=color_by_file_type[file_type],
                       label=file_type, alpha=0.7, width=width)

    # format
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    # remove x axis labels
    ax.set_xticks([])
    # log scale
    # ax.set_yscale('log')

    # custom legend for file ( types rectangles )
    handles = [plt.Rectangle((0, 0), 1, 1,
                             color=color_by_file_type[file_type],
                             alpha=0.7, label=file_type)
               for file_type in color_by_file_type.keys()]

    ax.legend(handles=handles, frameon=False, title='File Type')
    plt.tight_layout()

    # save
    plt.savefig(os.path.join(out, 'file_sizes.png'))

def main():
    args = parse_args()
    num_genes = 20386
    pval = 7.3

    tabix_times_file = plot_utils.read_tabix_times_data(args.tabix_files)
    new_times_file = plot_utils.read_new_times_data(args.new_files)
    tabix_times_genes, tabix_gene_hits = plot_utils.read_tabix_genes(args.tabix_genes)
    new_times_genes, new_gene_hits = plot_utils.read_new_genes(args.new_genes)


    file_names = list(tabix_times_file.keys())
    file_sizes = plot_utils.get_file_sizes(args.compressed_files_dir,
                                file_names)

    plot_compare_times_genes(tabix_times_genes,
                             new_times_genes,
                             tabix_gene_hits,
                             args.out)

    plot_compare_times_files(tabix_times_file,
                             new_times_file,
                             num_genes,
                             pval,
                             args.out)

    plot_file_sizes(file_sizes,
                    args.out)


if __name__ == '__main__':
    main()

