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
    parser.add_argument('--root', type=str, required=True,
                        help='dir with data')
    parser.add_argument('--bed', type=str, required=True,
                        help='bed file genes')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')
    return parser.parse_args()

def plot_compare_times_genes(file,
                             block_sizes,
                             codecs,
                             file_sizes,
                             data_dict_times, data_dict_records, data_dict_pval_hits,
                             gene_sizes,
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
                           figsize=(32, 18), dpi=300,
                           sharex=False, sharey=False)
    # fig.suptitle('Tabix vs. XXX\nsearch times by gene', fontweight='bold')
    for i, block_size in enumerate(block_sizes):
        for j, codec in enumerate(codecs):
            # x = tabix time, y = xxx time, color = number of hits
            x_data = []
            y_data = []
            pval_data = []
            hits_data = []
            gene_sizes_data = []
            instances_data = []
            avg_speedup = []
            indexes_data = []
            math_data = []
            tabix_wins = 0
            xxx_wins = 0
            for gene in data_dict_times[codec][block_size].keys():
                x_data.append(data_dict_times['tabix'][gene])
                y_data.append(data_dict_times[codec][block_size][gene])
                pval_data.append(data_dict_pval_hits[codec][block_size][gene])
                hits_data.append(data_dict_records[codec][block_size][gene])
                gene_sizes_data.append(gene_sizes[gene])
                # instances_data.append(gene_instances[gene])
                indexes_data.append(gene_indexes[gene])
                speedup = data_dict_times['tabix'][gene] / data_dict_times[codec][block_size][gene]
                avg_speedup.append(speedup)

                if data_dict_pval_hits[codec][block_size][gene] == 0:
                    math_data.append(0)
                else:
                    math_data.append(gene_indexes[gene])

                # count how many times tabix wins
                if data_dict_times['tabix'][gene] < data_dict_times[codec][block_size][gene]:
                    tabix_wins += 1
                elif data_dict_times['tabix'][gene] > data_dict_times[codec][block_size][gene]:
                    xxx_wins += 1
                else:
                    # tie
                    pass

            uncompressed_size = file_sizes['tsv']
            bgz_tbi_size = (file_sizes['bgz'] +
                               file_sizes['tbi'])
            xxx_size = (file_sizes[block_size][codec]['xxx'] +
                        file_sizes[block_size][codec]['gen'] +
                        file_sizes[block_size][codec]['pval'])

            # add text to plot which shows file size ratios
            comp_ratio = bgz_tbi_size / xxx_size
            # xxx_to_uncompressed = xxx_size / uncompressed_size
            text = 'Comp. Ratio: {:.4f}'.format(comp_ratio)
            ax[i, j].text(0.06, 0.95, text,
                          transform=ax[i, j].transAxes,
                          horizontalalignment='left',
                          verticalalignment='top',
                          fontsize=10,
                          bbox=dict(facecolor='none', alpha=0.5, edgecolor='none'))

            # add text to plot which shows avg speedup search time
            avg_speedup_file = np.mean(avg_speedup)
            text = 'Avg Speedup: {:.4f}'.format(avg_speedup_file)
            ax[i, j].text(0.06, 0.85, text,
                          transform=ax[i, j].transAxes,
                          horizontalalignment='left',
                          verticalalignment='top',
                          fontsize=10,
                          bbox=dict(facecolor='none', alpha=0.5, edgecolor='none'))

            # add text to plot which shows percent xxx wins
            text = '% STABIX Wins: {:.4f}'.format((xxx_wins / (xxx_wins + tabix_wins)) * 100)
            ax[i, j].text(0.06, 0.75, text,
                          transform=ax[i, j].transAxes,
                          horizontalalignment='left',
                          verticalalignment='top',
                          fontsize=10,
                          bbox=dict(facecolor='none', alpha=0.5, edgecolor='none'))

            hits = ax[i, j].scatter(x_data, y_data, s=20, alpha=0.1, color='darkblue')

            # for g in curious:
            #     # plot a large date point
            #     ax[i, j].scatter(data_dict_times['tabix'][g],
            #                      data_dict_times[codec][block_size][g],
            #                         c=curious_colors[g], s=500, alpha=0.7)

            # draw line y = x
            ax[i, j].plot([0, 0.05], [0, 0.05], color='black', linestyle='--')

            # ax[i, j].set_xlabel('Tabix Search Time (s)', fontsize=10)
            # ax[i, j].set_ylabel('STABIX Search Time (s)', fontsize=10)

            # format
            ax[i, j].spines['top'].set_visible(False)
            ax[i, j].spines['right'].set_visible(False)

            # add colorbar
            # cbar = plt.colorbar(hits, ax=ax[i, j])
            # cbar.set_label('Number of XXX Blocks')


    # add text at top of each row with codec name
    for i, codec in enumerate(codecs):
        title = 'Codec Comb.:\n{}'.format(codec)
        ax[0, i].set_title(title, fontsize=16, fontweight='bold')
        ax[2, i].set_xlabel('Tabix Search Time (s)', fontsize=10)
    # add text to left of each row with block size
    for i, block_size in enumerate(block_sizes):
        text = 'Block Size: {}'.format(block_size)
        ax[i, 0].text(-0.3, 0.5, 'Block Size: ' + block_size, fontsize=16, fontweight='bold',
                      rotation=90, verticalalignment='center', horizontalalignment='center',
                      transform=ax[i, 0].transAxes)
        ax[i, 0].set_ylabel('STABIX Search Time (s)', fontsize=10)


    # plt.tight_layout()
    plt.savefig(out_png)


def plot_compare_times_files(tabix_times,
                             xxx_times,
                             num_genes,
                             pval,
                             out):
    # tabix times are in seconds (python time.time())
    # xxx times are in microseconds (c++ chrono::high_resolution_clock::now())

    # scatter: x = tabix time, y = xxx time, color = number of hits
    fig, ax = plt.subplots(1, 1, figsize=(20, 8), dpi=300)

    for gwas_file, (genes, snps, tabix_time) in tabix_times.items():
        if gwas_file not in xxx_times:
            continue
        xxx_time = xxx_times[gwas_file]
        num_hits = genes
        ax.scatter(tabix_time, xxx_time, alpha=0.7)

    ax.set_xlabel('Tabix Search Time (s)')
    ax.set_ylabel('XXX Search Time (s)')
    ax.set_title('Tabix vs. XXX Search Times', fontweight='bold')

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
    plt.savefig(os.path.join(out, 'tabix_vs_xxx_search_times.png'))


def main():
    args = parse_args()
    file_names = ['continuous-103220-both_sexes']
    # codecs = ['bz2', 'xz', 'zlib', 'combo-xzb']
    # block_sizes = ['1000', '2000', '10000']
    codecs = ['bz2', 'deflate', 'xz', 'zlib', 'zstd', 'combo-fbb', 'combo-xbb', 'combo-xzb']
    block_sizes = ['1000', '2000', '5000', '10000', 'map']

    gene_sizes = plot_utils.get_gene_size(args.bed)
    gene_instances = plot_utils.get_gene_instances(args.bed)

    data_dir = args.root + 'data/'
    tabix_dir = args.root + 'tabix_output/'
    out_dir = args.out

    data_dict_times = defaultdict(dict)
    data_dict_records = defaultdict(dict)
    data_dict_pval_hits = defaultdict(dict)

    file_sizes = defaultdict(dict)

    gene_locations = plot_utils.read_bed_file(args.bed)
    # genomic_index = (plot_utils.read_genomic_index
    #                  ('data/UKBB/continuous-103220-both_sexes_output/zlib_5000.gen.idx'))
    # gene_indexes = plot_utils.get_bp_index_from_genomic_index(genomic_index,
    #                                                           gene_locations)

    for file in file_names:
        tabix_file = os.path.join(tabix_dir, file + '_tabix_output.txt')
        (tabix_gene_times,
         tabix_gene_records,
         tabix_pval_hits) = plot_utils.read_genes(tabix_file, False)
        data_dict_times['tabix'] = tabix_gene_times[file]
        data_dict_records['tabix'] = tabix_gene_records[file]
        data_dict_pval_hits['tabix'] = tabix_pval_hits[file]

        # get size of files
        tsv = data_dir + file + '.tsv'
        bgz = data_dir + file + '.tsv.bgz'
        tbi = data_dir + file + '.tsv.bgz.tbi'
        file_sizes[file] = {'tsv': os.path.getsize(tsv),
                            'bgz': os.path.getsize(bgz),
                            'tbi': os.path.getsize(tbi)}
        for block_size in block_sizes:
            for codec in codecs:
                file_name = file + '_' + block_size + '_' + codec
                file_path = os.path.join(data_dir, file_name + '/')
                if not os.path.exists(file_path):
                    print('File does not exist:', file_path)
                    continue

                # get size of files for xxx, gen, and pval
                xxx_file = file_path + file_name + '.grlz'
                gen_file = os.path.join(file_path, 'genomic.idx')
                pval_file = os.path.join(file_path, 'pval.idx')
                try:
                    file_sizes[file][block_size][codec] = {'xxx': os.path.getsize(xxx_file),
                                                       'gen': os.path.getsize(gen_file),
                                                       'pval': os.path.getsize(pval_file)}
                except KeyError:
                    try:
                        file_sizes[file][block_size] = {codec: {'xxx': os.path.getsize(xxx_file),
                                                                'gen': os.path.getsize(gen_file),
                                                                'pval': os.path.getsize(pval_file)}}
                    except KeyError:
                        file_sizes[file] = {block_size: {codec: {'xxx': os.path.getsize(xxx_file),
                                                                'gen': os.path.getsize(gen_file),
                                                                'pval': os.path.getsize(pval_file)}}}



                print('Reading file:', file_path)
                (xxx_gene_times,
                 xxx_gene_records,
                 xxx_pval_hits) = plot_utils.read_genes(file_path + file + '.query', True)

                plot_utils.check_gene_records(tabix_gene_records[file],
                                              xxx_gene_records[file])

                try:
                    data_dict_times[codec][block_size] = xxx_gene_times[file]
                    data_dict_records[codec][block_size] = xxx_gene_records[file]
                    data_dict_pval_hits[codec][block_size] = xxx_pval_hits[file]

                except KeyError:
                    data_dict_times[codec] = {block_size: xxx_gene_times[file]}
                    data_dict_records[codec] = {block_size: xxx_gene_records[file]}
                    data_dict_pval_hits[codec] = {block_size: xxx_pval_hits[file]}

        # plot compare times for each file
        out_png = os.path.join(out_dir, file + '_compare_times.png')
        plot_compare_times_genes(file,
                                 block_sizes, codecs,
                                 file_sizes[file],
                                 data_dict_times,
                                 data_dict_records,
                                 data_dict_pval_hits,
                                 gene_sizes, gene_instances,
                                 out_png)

    # plot_compare_times_files(tabix_times_file,
    #                          xxx_times_file,
    #                          num_genes,
    #                          pval,
    #                          args.out)
    #
    # plot_file_sizes(file_sizes,
    #                 args.out)


if __name__ == '__main__':
    main()

