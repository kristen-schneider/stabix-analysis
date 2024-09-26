import argparse

import matplotlib.pyplot as plt
import os
from collections import defaultdict

import plot_utils as plot_utils

def parse_args():
    parser = argparse.ArgumentParser(description='Plot New Results')
    parser.add_argument('--tabix_times', type=str, required=True,
                        help='tabix times')
    parser.add_argument('--results_dir', type=str, required=True,
                        help='new search results')
    parser.add_argument('--compressed_dir', type=str, required=True,
                        help='compressed files')
    parser.add_argument('--colors', type=str, required=True,
                        help='colors for each codec')
    parser.add_argument('--bed', type=str, required=True,
                        help='bed file')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')
    return parser.parse_args()

def plot_compare_new_times(file_sizes,
                           tabix_times,
                           files_gene_times,
                           block_sizes,
                           codecs,
                           gene_sizes,
                           colors,
                           out):

    x_range = range(0, len(block_sizes)+1)
    x_labels = [str(block_size) for block_size in block_sizes]
    x_labels.append('tabix')
    x_marks_dict = {x_labels[i]: i for i in x_range}

    # offset each codec by 0.2 from x-mark
    codec_offsets = {'bz2': -0.2, 'xz': 0.0, 'zlib': 0.2}

    # x-axis is the block size, y-axis is the time

    tabix_plotted = False
    fig, ax = plt.subplots(2, 1, figsize=(10, 12), dpi=300)
    for file, gene_times in files_gene_times.items():
        codec, block_size, new = file.split('_')

        x = x_marks_dict[block_size] + codec_offsets[codec]
        for trait, genes in gene_times.items():
            # for gene, time in genes.items():
            #     if time > .1:
            #         print(codec, block_size, gene, time, gene_sizes[gene])

            # plot the tabix time only once
            if not tabix_plotted:
                tabix_times_curr = tabix_times[trait]
                tabix_vp = ax[0].violinplot(tabix_times_curr.values(),
                                            positions=[x_marks_dict['tabix']],
                                            showmeans=False,
                                            showmedians=True,
                                            widths=0.2)
                # color the violins fill and lines around them with the codec color
                plt.setp(tabix_vp['bodies'], facecolor=colors['tabix'], edgecolor=colors['tabix'])
                plt.setp(tabix_vp['cmedians'], color=colors['tabix'])
                plt.setp(tabix_vp['cbars'], color=colors['tabix'])
                plt.setp(tabix_vp['cmins'], color=colors['tabix'])
                plt.setp(tabix_vp['cmaxes'], color=colors['tabix'])
                tabix_plotted = True

                # write the gene name with the largest time
                max_time = max(tabix_times_curr.values())
                max_gene = [gene for gene, time in tabix_times_curr.items() if time == max_time][0]
                ax[0].text(x_marks_dict['tabix'], max_time, max_gene, fontsize=6)
                # write mean search time
                mean_time = sum(tabix_times_curr.values()) / len(tabix_times_curr)
                ax[0].text(x_marks_dict['tabix'], mean_time, 'Mean: {:.4f}'.format(mean_time), fontsize=6)


                # plot file size for tabix
                file_size_name = trait + '_' + block_size + '_' + codec
                file_size = file_sizes[file_size_name]
                ax[1].bar(x_marks_dict['tabix'], file_size['bgz'] + file_size['tbi'],
                          color=colors['tabix'], width=0.1)

            new_vp = ax[0].violinplot(genes.values(),
                                      positions=[x],
                                      showmeans=False,
                                      showmedians=True,
                                      widths=0.2)

            # write the gene name with the largest time
            max_time = max(genes.values())
            max_gene = [gene for gene, time in genes.items() if time == max_time][0]
            ax[0].text(x, max_time, max_gene, fontsize=6)
            # write mean search time
            mean_time = sum(genes.values()) / len(genes)
            ax[0].text(x, mean_time, 'Mean: {:.4f}'.format(mean_time), fontsize=6)


            # color the violins fill and lines around them with the codec color
            plt.setp(new_vp['bodies'], facecolor=colors[codec], edgecolor=colors[codec])
            plt.setp(new_vp['cmedians'], color=colors[codec])
            plt.setp(new_vp['cbars'], color=colors[codec])
            plt.setp(new_vp['cmins'], color=colors[codec])
            plt.setp(new_vp['cmaxes'], color=colors[codec])


            # plot file size on second plot
            file_size_name = trait + '_' + block_size + '_' + codec
            file_size = file_sizes[file_size_name]
            ax[1].bar(x, file_size['new'] + file_size['gen'] + file_size['pval'],
                      color=colors[codec], width=0.1)

    # rectagle legend for codecs
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[codec], label=codec) for codec in codecs]
    # add black one for tabix
    handles.append(plt.Rectangle((0, 0), 1, 1, color=colors['tabix'], label='tabix'))
    ax[0].legend(handles=handles, title='Codec', loc='upper right', frameon=False)

    ## add text box with number of genes
    num_genes = 20385
    # no border on text box
    ax[0].text(0.2, 0.7, 'Number of Genes Queried:\n {}'.format(num_genes),
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax[0].transAxes,
            bbox=dict(facecolor='none', edgecolor='none'))
    ax[0].set_xticks(x_range)
    ax[0].set_xticklabels(x_labels)
    ax[0].set_xlabel('Block Size')
    ax[0].set_ylabel('Decompression Time (s)')
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)

    ax[1].set_xticks(x_range)
    ax[1].set_xticklabels(x_labels)
    ax[1].set_xlabel('Block Size')
    ax[1].set_ylabel('File Size (GB)')
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)



    plt.tight_layout()
    plt.savefig(os.path.join(out, 'new_search_times_investigation.png'))



def main():
    args = parse_args()
    results_dir = args.results_dir
    out = args.out

    codecs = ['bz2', 'xz']
    block_sizes = [5000, 10000]
    file_names = ['continuous-103220-both_sexes']
    suffix = '_new.txt'
    ext = '.txt'

    colors = plot_utils.read_colors(args.colors)
    gene_sizes = plot_utils.get_gene_size(args.bed)

    # sort by size of gene
    # sorted_gene_sizes = sorted(gene_sizes.items(), key=lambda x: x[1], reverse=True)
    # print largest genes
    # print(sorted_gene_sizes[:10])

    all_files = []
    for codec in codecs:
        for block_size in block_sizes:
            all_files.append(os.path.join(results_dir, '{}_{}'.format(codec, block_size) + suffix))

    files_gene_times = {}
    file_gene_records = {}
    for file in all_files:
        file_basename = os.path.basename(file).replace(ext, '')
        files_gene_times[file_basename], file_gene_records[file_basename] = plot_utils.read_new_genes(file)

    tabix_times, tabix_hits = plot_utils.read_tabix_genes(args.tabix_times)

    file_sizes = plot_utils.get_file_sizes(args.compressed_dir,
                                           file_names,
                                           codecs,
                                           block_sizes)

    plot_compare_new_times(file_sizes,
                           tabix_times,
                           files_gene_times,
                           block_sizes,
                           codecs,
                           gene_sizes,
                           colors,
                           out)




if __name__ == '__main__':
    main()