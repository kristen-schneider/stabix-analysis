import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from pandas.compat.numpy.function import validate_take_with_convert
from pysam.libcvcf import defaultdict

import plot_utils as plot_utils


def parse_args():
    parser = argparse.ArgumentParser(description='Plot Tabix vs. XXX Results')
    parser.add_argument('--data', type=str, required=True,
                        help='dir with data')
    parser.add_argument('--bed', type=str, required=True,
                        help='bed file genes')
    parser.add_argument('--comp_dir', type=str, required=True,
                        help='dir with compressed files')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')
    parser.add_argument('--compressed_files_dir', type=str, required=True,
                        help='directory with compressed files')
    return parser.parse_args()

def plot_compare_times_genes(block_sizes,
                             codecs,
                             file_sizes,
                             data_dict_times, data_dict_records, data_dict_pval_hits,
                             gene_sizes,
                             gene_instances,
                             gene_indexes,
                             out_png):

    colormap = 'viridis'
    colormap_bins = [0, 1, 2, 3, 4, 5, 20, 60]
    cmap = plt.get_cmap(colormap)
    norm = mcolors.BoundaryNorm(colormap_bins, cmap.N)
    curious = ['MAN2A2', 'H1-0', 'SLC39A9']
    curious_colors = {'MAN2A2': 'green',
                      'H1-0': 'yellow',
                      'SLC39A9': 'red'}

    fig, ax = plt.subplots(len(block_sizes), len(codecs),
                           figsize=(25, 20), dpi=300,
                           sharex=False, sharey=False)
    # fig.suptitle('Tabix vs. XXX\nsearch times by gene', fontweight='bold')
    for i, block_size in enumerate(block_sizes):
        for j, codec in enumerate(codecs):
            # x = tabix time, y = lzr time, color = number of hits
            x_data = []
            y_data = []
            pval_data = []
            hits_data = []
            gene_sizes_data = []
            instances_data = []
            indexes_data = []
            math_data = []
            tabix_wins = 0
            lzr_wins = 0
            for gene in data_dict_times[codec][block_size].keys():
                x_data.append(data_dict_times['tabix'][gene])
                y_data.append(data_dict_times[codec][block_size][gene])
                pval_data.append(data_dict_pval_hits[codec][block_size][gene])
                hits_data.append(data_dict_records[codec][block_size][gene])
                gene_sizes_data.append(gene_sizes[gene])
                instances_data.append(gene_instances[gene])
                indexes_data.append(gene_indexes[gene])

                if data_dict_pval_hits[codec][block_size][gene] == 0:
                    math_data.append(0)
                else:
                    math_data.append(gene_indexes[gene])

                # count how many times tabix wins
                if data_dict_times['tabix'][gene] < data_dict_times[codec][block_size][gene]:
                    tabix_wins += 1
                elif data_dict_times['tabix'][gene] > data_dict_times[codec][block_size][gene]:
                    lzr_wins += 1
                else:
                    # tie
                    pass

                # if (data_dict_times['tabix'][gene] < 0.06
                #         and data_dict_times[codec][block_size][gene] > 0.035
                #         and block_size == '2000' and codec == 'bz2'
                #         and (data_dict_pval_hits[codec][block_size][gene] * gene_indexes[gene]) <= 20):
                #     print('Gene:', gene)
                #     # print('-Hits:', data_dict_records[codec][block_size][gene])
                #     print('-P-value Hits:', data_dict_pval_hits[codec][block_size][gene])
                #     # print('-Records:', gene_instances[gene])
                #     print('-Blocks:', gene_indexes[gene])


                # # print genes where lzr time is greater than tabix time
                # # for best block size
                # if (data_dict_times['tabix'][gene] < 0.025
                #         and data_dict_times[codec][block_size][gene] > 0.035
                #         and block_size == '2000'
                #         and codec == 'bz2'
                #         and data_dict_pval_hits[codec][block_size][gene] == 0):
                #     print('Gene:', gene)
                #     print('-Tabix Time:', data_dict_times['tabix'][gene])
                #     print('-LZR Time:', data_dict_times[codec][block_size][gene])
                #     print('-Gene Hits:', data_dict_records[codec][block_size][gene])
                #     print('-P-value Hits:', data_dict_pval_hits[codec][block_size][gene])

            uncompressed_size = file_sizes[block_size][codec]['tsv']
            compressed_size = (file_sizes[block_size][codec]['bgz'] +
                               file_sizes[block_size][codec]['tbi'])
            lzr_size = (file_sizes[block_size][codec]['new'] +
                        file_sizes[block_size][codec]['gen'] +
                        file_sizes[block_size][codec]['pval'])

            # add text to plot which shows file size ratios
            lzr_to_tabix = lzr_size / compressed_size
            # lzr_to_uncompressed = lzr_size / uncompressed_size
            text = 'XXX/Tabix (comp ratio):\n{:.4f}'.format(lzr_to_tabix)
            ax[i, j].text(0.25, 0.95, text,
                          transform=ax[i, j].transAxes,
                          horizontalalignment='center',
                          verticalalignment='top',
                          fontsize=10,
                          bbox=dict(facecolor='none', alpha=0.5, edgecolor='black'))

            # add text to plot which shows how many times tabix wins
            text = 'Tabix Wins: {}\nXXX Wins: {}'.format(tabix_wins, lzr_wins)
            ax[i, j].text(0.75, 0.95, text,
                            transform=ax[i, j].transAxes,
                            horizontalalignment='center',
                            verticalalignment='top',
                            fontsize=10,
                            bbox=dict(facecolor='none', alpha=0.5, edgecolor='black'))

            hits = ax[i, j].scatter(x_data, y_data, s=20, alpha=0.7)

            for g in curious:
                # plot a large date point
                ax[i, j].scatter(data_dict_times['tabix'][g],
                                 data_dict_times[codec][block_size][g],
                                    c=curious_colors[g], s=500, alpha=0.7)

            # draw line y = x
            ax[i, j].plot([0, 0.1], [0, 0.1], color='black', linestyle='--')

            ax[i, j].set_xlabel('Tabix Search Time (s)')
            ax[i, j].set_ylabel('LZR Search Time (s)')

            # format
            ax[i, j].spines['top'].set_visible(False)
            ax[i, j].spines['right'].set_visible(False)
            # set x and y limits
            ax[i, j].set_xlim(0, 0.18)
            ax[i, j].set_ylim(0, 0.18)
            # add colorbar
            # cbar = plt.colorbar(hits, ax=ax[i, j])
            # cbar.set_label('Number of XXX Blocks')


    # add text at top of each row with codec name
    for i, codec in enumerate(codecs):
        title = 'Codec: {}'.format(codec)
        ax[0, i].set_title(title, fontsize=10, fontweight='bold')
    # add text to left of each row with block size
    for i, block_size in enumerate(block_sizes):
        text = 'Block Size: {}'.format(block_size)
        ax[i, 0].text(-0.3, 0.5, 'Block Size: ' + block_size, fontsize=16, fontweight='bold',
                      rotation=90, verticalalignment='center', horizontalalignment='center',
                      transform=ax[i, 0].transAxes)

    # plt.tight_layout()
    plt.savefig(out_png)


def plot_compare_times_files(tabix_times,
                             new_times,
                             num_genes,
                             pval,
                             out):
    # tabix times are in seconds (python time.time())
    # new times are in microseconds (c++ chrono::high_resolution_clock::now())

    # scatter: x = tabix time, y = new time, color = number of hits
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), dpi=300)

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


def main():
    args = parse_args()
    num_genes = 20386
    pval = 7.3

    # file_names = ['continuous-103220-both_sexes']
    file_names = ['continuous-103220-both_sexes', 'continuous-30130-both_sexes-irnt']
    codecs = ['bz2', 'xz', 'combo']
    block_sizes = ['2000', '5000', '10000']
    gene_sizes = plot_utils.get_gene_size(args.bed)
    gene_instances = plot_utils.get_gene_instances(args.bed)

    data_dir = args.data
    out_dir = args.out

    data_dict_times = defaultdict(dict)
    data_dict_records = defaultdict(dict)
    data_dict_pval_hits = defaultdict(dict)

    gene_locations = plot_utils.read_bed_file(args.bed)
    genomic_index = (plot_utils.read_genomic_index
                     ('data/UKBB/continuous-103220-both_sexes_output/zlib_5000.gen.idx'))
    gene_indexes = plot_utils.get_bp_index_from_genomic_index(genomic_index,
                                                              gene_locations)

    file_sizes = plot_utils.get_file_sizes(args.comp_dir,
                                           file_names, codecs,
                                           block_sizes)

    for file in file_names:
        file_dir = os.path.join(data_dir, file + '_output')
        tabix_file = os.path.join(file_dir, 'tabix_output.txt')
        tabix_gene_times, tabix_gene_records, tabix_pval_hits = plot_utils.read_genes(tabix_file, False)
        data_dict_times['tabix'] = tabix_gene_times[file]
        data_dict_records['tabix'] = tabix_gene_records[file]
        data_dict_pval_hits['tabix'] = tabix_pval_hits[file]
        for codec in codecs:
            for block_size in block_sizes:
                file_name = codec + '_' + str(block_size) + '.txt'
                file_path = os.path.join(file_dir, file_name)

                if not os.path.exists(file_path):
                    print('File does not exist:', file_path)
                    continue

                print('Reading file:', file_path)
                lzr_gene_times, lzr_gene_records, lzr_pval_hits = plot_utils.read_genes(file_path, True)
                plot_utils.check_gene_records(tabix_gene_records[file],
                                              lzr_gene_records[file])

                try:
                    data_dict_times[codec][block_size] = lzr_gene_times[file]
                    data_dict_records[codec][block_size] = lzr_gene_records[file]
                    data_dict_pval_hits[codec][block_size] = lzr_pval_hits[file]

                except KeyError:
                    data_dict_times[codec] = {block_size: lzr_gene_times[file]}
                    data_dict_records[codec] = {block_size: lzr_gene_records[file]}
                    data_dict_pval_hits[codec] = {block_size: lzr_pval_hits[file]}

        # plot compare times for each file
        out_png = os.path.join(out_dir, file + '_compare_times.png')
        plot_compare_times_genes(block_sizes, codecs,
                                 file_sizes[file],
                                 data_dict_times,
                                 data_dict_records,
                                 data_dict_pval_hits,
                                 gene_sizes, gene_instances, gene_indexes,
                                 out_png)

    # plot_compare_times_files(tabix_times_file,
    #                          new_times_file,
    #                          num_genes,
    #                          pval,
    #                          args.out)
    #
    # plot_file_sizes(file_sizes,
    #                 args.out)


if __name__ == '__main__':
    main()

