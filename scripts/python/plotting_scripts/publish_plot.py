from cProfile import label
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
import pandas as pd
from PIL.ImageColor import colormap
from scipy.constants import alpha
from setuptools.command.rotate import rotate
from matplotlib import cm
from matplotlib import colors as mcolors

import plot_utils as pltut

def plot_scatter_sizes(file_sizes,
                       out):
    # for each file, plot:
    # tsv (x-axis), bgz +tbi (y-axis),
    # tsv (x-axis), xxx + gen + pval (y-axis)
    # in a scatter plot

    fig, ax = plt.subplots(1, 1, figsize=(10, 5), dpi=300)
    colors = {'bgz': 'olivedrab', 'tbi': 'green',
              'xxx': 'darkorange', 'gen': 'orange'}

    tabix_x = []
    bgz_y = []
    tbi_y = []
    xxx_y = []
    xxxidx_y = []

    for f in file_sizes:
        tabix_x.append(file_sizes[f]['tsv']/1e9)
        bgz_y.append((file_sizes[f]['bgz'] + file_sizes[f]['tbi'])/1e9)
        # tbi_y.append((file_sizes[f]['tbi']/1e9))
        xxx_y.append((file_sizes[f]['xxx'] + file_sizes[f]['gen'] + file_sizes[f]['pval'])/1e9)
        # xxxidx_y.append((file_sizes[f]['gen']/1e9 + file_sizes[f]['pval']/1e9))

    ax.plot(tabix_x, bgz_y, color=colors['bgz'], label='bgz', marker='o', alpha=0.5)
    # ax.scatter(tabix_x, tbi_y, color=colors['tbi'], label='tbi')
    ax.plot(tabix_x, xxx_y, color=colors['xxx'], label='xxx', marker='o', alpha=0.5)
    # ax.scatter(tabix_x, xxxidx_y, color=colors['gen'], label='gidx + pidx')

    ax.set_xlabel('Uncompressed TSV File Size (GB)')
    ax.set_ylabel('Compressed File Size (GB)')

    # custom legend
    # a line for each color in colors
    # a label for each color in colors
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['bgz'], markersize=10, label='bgz'),
                          # plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['tbi'], markersize=10, label='tbi'),
                            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['xxx'], markersize=10, label='xxx')]
                            # plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['gen'], markersize=10, label='gen + pval')]
    ax.legend(legend_elements, ['bgz + tbi', 'STABIX + gidx + sidx'], loc='upper left', title='Compression Types', frameon=False)

    # spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(out + 'file_sizes_scatter.png')


def plot_bar_sizes(file_sizes,
                   out):
    # for each file, plot:
    # tsv (left), bgz +tbi (middle), xxx + gen + pval (right)
    # in a bar plot
    # x-axis: file names
    # y-axis: file sizes
    # colors: tsv = blue, bgz = green, tbi = red, xxx = purple, gen = orange, pval = brown
    # legend: tsv, bgz, tbi, xxx, gen, pval
    # title: File sizes
    # save to out + 'file_sizes.png'

    fig, ax = plt.subplots(1, 1, figsize=(8, 3.5), dpi=300)

    x_centers = list(range(len(file_sizes)))
    x_labels = file_sizes.keys()

    file_name_traits = {'continuous-103220-both_sexes': 'shellfish\nintake',
                        'phecode-282.5-both_sexes': 'sickle cell\nanemia',
                        'categorical-20096-both_sexes-2': 'size of\nred wine\nglass drunk',
                        'phecode-696.4-both_sexes': 'psoriasis',
                        'continuous-50-both_sexes-irnt': 'standing\nheight',
                        'continuous-30130-both_sexes-irnt': 'monocyte\ncount',
                        'phecode-250.2-both_sexes': 'type 2\ndiabetes',
                        'phecode-250-both_sexes': 'diabetes\nmellitus',
                        'categorical-1210-both_sexes-1210': 'snoring',
                        'categorical-20116-both_sexes-0': 'smoking\nstatus'}

    colors = {'tsv': 'grey',
              'bgz': 'darkorange', 'tbi': 'red',
              'xxx': 'darkblue', 'gen': 'red', 'pval': 'red'}

    for f in file_sizes:
        x = x_centers.pop(0)
        x_tsv = x - 0.2
        x_bgz = x
        x_xxx = x + 0.2
        annotation_pad = 0.3

        tsv = ax.bar(x_tsv, file_sizes[f]['tsv']/1e9, color=colors['tsv'], width=0.2, label='tsv')
        tabix = ax.bar(x_bgz, file_sizes[f]['bgz']/1e9 + file_sizes[f]['tbi']/1e9, color=colors['bgz'], width=0.2, label='bgz + tbi')
        xxx = ax.bar(x_xxx, file_sizes[f]['xxx']/1e9 + file_sizes[f]['gen']/1e9 + file_sizes[f]['pval']/1e9, color=colors['xxx'], width=0.2, label='xxx + gen + pval')

        # annotate the bars with file size rotated 90 degrees and add padding from the top of the bar
        ax.text(x_tsv, file_sizes[f]['tsv']/1e9 + annotation_pad, '{:.3f}'.format(file_sizes[f]['tsv']/1e9),
                ha='center', fontsize=8, rotation=90)
        ax.text(x_bgz, file_sizes[f]['bgz']/1e9 + file_sizes[f]['tbi']/1e9 + annotation_pad, '{:.3f}'.format(file_sizes[f]['bgz']/1e9 + file_sizes[f]['tbi']/1e9),
                ha='center', fontsize=8, rotation=90)
        ax.text(x_xxx, file_sizes[f]['xxx']/1e9 + file_sizes[f]['gen']/1e9 + file_sizes[f]['pval']/1e9 + annotation_pad, '{:.3f}'.format(file_sizes[f]['xxx']/1e9 + file_sizes[f]['gen']/1e9 + file_sizes[f]['pval']/1e9),
                ha='center', fontsize=8, rotation=90)

        # # add compression ratio annotation
        # compression_ratio = ((file_sizes[f]['bgz']/1e9 + file_sizes[f]['tbi']/1e9) /
        #                      (file_sizes[f]['xxx']/1e9 + file_sizes[f]['gen']/1e9 + file_sizes[f]['pval']/1e9))
        # ax.text(x_bgz + 0.2, (file_sizes[f]['bgz']/1e9 + file_sizes[f]['tbi']/1e9) + annotation_pad + 2, '{:.3f}'.format(compression_ratio),
        #         ha='center', fontsize=8, fontweight='bold')


    ax.set_xticks(list(range(len(file_sizes))))
    ax.set_xticklabels([file_name_traits[f] for f in file_sizes], fontsize=7)
    ax.set_xlabel('File Traits')
    ax.set_ylabel('File Size (GB)')

    # custom legend
    # a rectangle for each color in colors
    # a label for each color in colors
    legend_elements = [plt.Rectangle((0, 0), 1, 1, fc=colors['tsv']),
                          plt.Rectangle((0, 0), 1, 1, fc=colors['bgz']),
                          plt.Rectangle((0, 0), 1, 1, fc=colors['xxx'])]
    labels = ['tsv (uncompressed)', f'bgzip + {name}', 'STABIX + gidx + sidx']
    ax.legend(legend_elements, labels, loc='upper left',
              title='Compression Types',
              fontsize=7,
              frameon=False)


    fig.suptitle('Compressed File Sizes', fontsize=12, fontweight='bold')

    # spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(out + 'file_sizes_bar.png')

def plot_all_files_gene_times(all_gene_times,
                                all_gene_records,
                                all_gene_pval_hits,
                                all_genomic_indexes,
                                file_sizes,
                                out_names,
                                block_sizes,
                                genes,
                                out,
                                name='tabix'):

    option = 'snps'

    # create main axies, leaving 25% of the space at top and on right for the histograms
    # fig, ax1 = plt.subplots(1, 1,
    #                         figsize=(8, 5), dpi=300)
    fig = plt.figure(figsize=(8, 8), dpi=300)
    # ax1 = fig.add_gridspec(1, 1, left=0.1, right=0.9, top=0.9, bottom=0.1).subplots()
    ax1 = fig.add_gridspec(1, 1, left=0.1, right=.75, top=.75, bottom=0.1).subplots()
    # ax1.set(aspect=0.55)


    # for each file, plot:
    # x-axis: tabix gene times
    # y-axis: gene names
    # file_size_data = []
    custom_colorbar = 'plasma'
    cmap = plt.get_cmap(custom_colorbar)

    sizes_bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    blocks_bins = [0, 1, 2, 3, 4, 5, 10, 20, 25, 30]
    pval_bins = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000]


    x_data = []
    y_data = []
    file_sizes_data = []
    pval_hits_data = []
    records_data = []
    num_blocks_data = []
    math_data = [] # (pval_hits_data / num_blocks_data) * records_data

    plot_options = {'none': {'save': 'gene_times.png', 'bins': None, 'colorbar': None},
                    'size': {'save': 'file_sizes.png', 'bins': sizes_bins, 'colorbar': 'Uncompressed File Size (GB)'},
                    'blocks': {'save': 'num_blocks.png', 'bins': blocks_bins,
                               'colorbar': 'Num XXX Blocks To Decompress'},
                    'pval': {'save': 'pval_hits.png', 'bins': pval_bins, 'colorbar': 'Num pval hits'},
                    'snps': {'save': 'snps_hits.png', 'bins': pval_bins, 'colorbar': 'Num sig SNPs'}}

    for f in all_gene_times:
        x = all_gene_times[f]['tabix']
        for block_size in block_sizes:
            for name in out_names:
                file_size = -1
                # if (f == 'phecode-282.5-both_sexes'):
                #         # or f == 'icd10-E11-both_sexes'
                #         or f =='phecode-282.5-both_sexes'):
                # if 'continuous-50-both_sexes-irnt' not in f: #!= 'icd10-E11-both_sexes':
                #     continue
                print(f)

                try:
                    curr_genomic_index = all_genomic_indexes[f][block_size][name]
                    curr_bp_indexes = pltut.get_bp_index_from_genomic_index(curr_genomic_index, genes)
                except KeyError:
                    continue
                try:
                    y = all_gene_times[f][block_size][name]
                    for g in genes:
                        x_data.append(x[g])
                        y_data.append(y[g])
                        file_size = file_sizes[f]['tsv'] / 1e9
                        file_sizes_data.append(file_size)
                        pval_hits_data.append(all_gene_pval_hits[f][block_size][name][g])
                        records_data.append(all_gene_records[f][block_size][name][g])
                        num_blocks_data.append(curr_bp_indexes[g])
                        # (records_data / num_blocks_data)

                        plot_options_gene = {'none': g,
                                             'size': g,
                                             'blocks': curr_bp_indexes[g],
                                             'pval': all_gene_pval_hits[f][block_size][name][g],
                                             'snps': all_gene_records[f][block_size][name][g]}

                        # if y[g] > 0.1 and y[g] > x[g]:
                        #     # label gene
                        #     ax1.text(x[g], y[g]+0.005, plot_options_gene[option],
                        #              fontsize=6, ha='center', va='center')

                    #     # label = file_sizes[f]['tsv'] / 1e9 to 3 decimal places
                    #     # label = '{:.3f}'.format(file_sizes[f]['tsv'] / 1e9)

                except KeyError:
                    continue

    plot_options_data = {'size': file_sizes_data,
                         'blocks': num_blocks_data,
                         'pval': pval_hits_data,
                         'snps': records_data}

    print('done getting data...')
    print('plotting scatter plot...')

    # plot scatter plot
    # hb = ax1.hexbin(x_data, y_data, gridsize=50, cmap=cmap,
    #                mincnt=1, bins='log')

    sct = ax1.scatter(x_data, y_data,
                      label=label, alpha=0.1,
                      # color='cadetblue')
                      ## c = plot_options_data[option], cmap=custom_colorbar)
                      c=plot_options_data[option],
                      cmap=custom_colorbar,
                      norm=mcolors.BoundaryNorm(plot_options[option]['bins'], cmap.N, clip=True))

    # plot histogram of xxx times on right rotated
    # left=0.1, right=0.75, top=.9, bottom=0.1
    # 0.75: This is the left position of the inset axes, meaning the inset starts at 75% of the width of the figure from the left.
    # 0: This is the bottom position of the inset axes, meaning it starts at the bottom of the figure (0%).
    # 0.25: This is the width of the inset axes, meaning it occupies 25% of the figure's width.
    # 1: This is the height of the inset axes, meaning it occupies the full height of the figure (100%).
    print('plotting histograms...')
    ax_histy = ax1.inset_axes([1, 0, .25, 1], sharey=ax1)
    ax_histy.hist(y_data, bins=100,
                  orientation='horizontal',
                  color='darkblue')
    ax_histy.set_xticks([])
    plt.setp(ax_histy.get_yticklabels(), visible=False)
    ax_histy.set_xscale('log')
    ax_histy.set_xticks([10, 1000, 100000])
    ax_histy.spines['top'].set_visible(False)
    ax_histy.spines['right'].set_visible(False)

    # plot histograms of tabix times on top
    ax_histx = ax1.inset_axes([0, 1.0, 1, 0.25], sharex=ax1)
    ax_histx.hist(x_data, bins=100,
                    color='darkblue')
    ax_histx.set_yticks([])
    plt.setp(ax_histx.get_xticklabels(), visible=False)
    ax_histx.set_yscale('log')
    ax_histx.set_yticks([10, 1000, 100000])
    ax_histx.spines['top'].set_visible(False)
    ax_histx.spines['right'].set_visible(False)

    print('plotting colorbar...')
    # plot colorbar, scale to the right of the scatter plot
    # make colorbar
    # cbar = fig.colorbar(sct, ax=ax_histy)
    cbar = fig.colorbar(ScalarMappable(cmap=sct.get_cmap(), norm=sct.norm), ax=ax_histy)
    # cbar = fig.colorbar(sct)
    #                     # orientation='vertical')#,
    #                     # pad=-0.15,
    #                     # shrink=0.8)
    #                     # label='counts in bin')
    # # cbar.set_label('counts in bin', fontsize=12)
    cbar.set_label(plot_options[option]['colorbar'], fontsize=8)  # Change the font size here
    # # cbar.ax.tick_params(labelsize=8)
    cbar.ax.yaxis.set_label_coords(-1., 0.5)
    # # ## left, bottom, width, height
    cbar.ax.set_position([0.65, 0.25, 0.3, 0.4])
    # cbar.ax.set_alpha(0.9)

    # plot a line with slope 1
    ax1.plot([0, 0.12], [0, 0.12], color='black', linestyle='--')

    # make x and y scale the same
    ax1.set_xlim([0, 0.24])
    ax1.set_ylim([0, 0.24])

    # add textbox which shows the block size and name
    string_label = ('Block size: ' + block_size + '\n'
                    + 'Comp. Type: ' + name + '\n'
                    + 'Num Files: ' + str(10) + '\n'
                    + 'Num Genes: ' + str(len(genes)))
    ax1.text(0.03, 0.9, string_label,
            transform=ax1.transAxes, fontsize=10,
            verticalalignment='top')

    # # add textbox with wins data
    # ax1.text(0.4, 0.9, 'File Size(GB), %XXX wins',
    #         transform=ax1.transAxes, fontsize=10, fontweight='bold',
    #         verticalalignment='top')
    # # ax1.text(0.4, 0.9, 'Total sig SNPs, %XXX wins',
    # #         transform=ax1.transAxes, fontsize=10, fontweight='bold',
    # #         verticalalignment='top')
    # y_pos = 0.85
    # for f in wins:
    #     ax1.text(0.5, y_pos, str('{:.3f}'.format(file_sizes[f]['tsv']/1e9)) + ', ' + str('{:.3f}'.format(wins[f])),
    #             transform=ax1.transAxes, fontsize=10,
    #             verticalalignment='top')
    #     # avg_num_records = sum(all_gene_records[f]['2000']['combo'].values())# / len(all_gene_records[f]['2000']['combo'])
    #     # ax1.text(0.5, y_pos, str('{:.0f}'.format(avg_num_records)) + ', ' + str('{:.3f}'.format(wins[f])),
    #     #         transform=ax1.transAxes, fontsize=10,
    #     #         verticalalignment='top')
    #     y_pos -= 0.05

    ax1.set_xlabel(f'{name.capitalize()} decompression time (s)')
    ax1.set_ylabel('STABIX decompression time (s)')
    # fig.suptitle('Gene Decompression Times', fontsize=12, fontweight='bold')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # plt.tight_layout()
    plt.savefig(out + plot_options[option]['save'])
    # plt.savefig(out + 'hexbin-' + plot_options[option]['save'])
    # plt.savefig(out + 'height_' + plot_options[option]['save'])



