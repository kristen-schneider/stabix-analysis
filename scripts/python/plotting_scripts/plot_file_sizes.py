from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.pyplot import annotate
from setuptools.command.rotate import rotate


def plot_file_sizes_scatter(gzip_sizes,
                    bgzip_sizes,
                    kzip_sizes,
                    codec_colors,
                    out_png):
    # get average gzip size
    avg_gzip_size = np.mean(gzip_sizes)
    # get average bgzip size
    avg_bgzip_size = np.mean(bgzip_sizes)

    # get average kzip sizes for each block size and codec
    avg_kzip_sizes = {}
    for block_size in kzip_sizes:
        avg_kzip_sizes[block_size] = {}
        for codec in kzip_sizes[block_size]:
            avg_kzip_sizes[block_size][codec] = np.mean(kzip_sizes[block_size][codec])

    # plot the average sizes
    fig, ax = plt.subplots(1, 1, figsize=(15, 7), dpi=300)
    # y = size in bytes
    # x = block size
    # color = codec cocktail

    # plot the average gzip size in red horizontal line
    plt.axhline(y=avg_gzip_size,
                color='black', linestyle='solid', linewidth=3,
                label='gzip')

    # plot the average bgzip size in blue horizontal line
    plt.axhline(y=avg_bgzip_size,
                color='black', linestyle='dashed', linewidth=3,
                label='bgzip')

    # single point for kzip codec cocktail sizes
    x_data = range(len(kzip_sizes.keys()))
    x_idx = 0
    for block_size in avg_kzip_sizes:
        x_point = x_data[x_idx]
        for codec in avg_kzip_sizes[block_size]:
            # jitter the x values
            x_point += np.random.uniform(-0.15, 0.15)
            plt.scatter(x_point,
                        avg_kzip_sizes[block_size][codec],
                        color=codec_colors[codec],
                        s=200,
                        marker='o',
                        edgecolor='black',
                        alpha=0.5)

        x_idx += 1

    # title
    plt.title('Compressed File Size\nby Block Size and Codec',
              fontsize=16, fontweight='bold')


    # FORMATTING
    # log scale
    plt.yscale('log')
    # legend, no frame
    # add custom legend for codecs
    codec_legend_elements = [Patch(facecolor=codec_colors[codec], edgecolor='black', label=codec)
                            for codec in codec_colors]
    # add bgzip dashed line to legend
    codec_legend_elements.insert(0, Line2D([0], [0],
                                           color='black', linestyle='dashed', linewidth=3,
                                           label='bgzip'))
    # add gzip solid line to legend
    codec_legend_elements.insert(1, Line2D([0], [0],
                                           color='black', linestyle='solid', linewidth=3,
                                           label='gzip'))

    plt.legend(handles=codec_legend_elements,
               title='Compression Method', loc='best', frameon=False,
               fontsize=12, title_fontsize=14)

    # label the y-axis
    plt.ylabel('Compressed file size\n(bytes)', fontsize=14)

    # label the x-axis
    x_ticks = range(len(kzip_sizes.keys()))
    x_tick_labels = [f'{block_size}' for block_size in kzip_sizes.keys()]
    plt.xticks(x_ticks, x_tick_labels, fontsize=12)
    plt.xlabel('Block size\n(number of records in block)', fontsize=14)

    # remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # SAVE
    plt.savefig(out_png)


def plot_file_sizes_violin(gzip_sizes,
                        bgzip_sizes,
                        kzip_sizes,
                        kzip_genomic_idx_sizes,
                        kzip_pvalue_idx_sizes,
                        codec_colors,
                        out_png):

    scale = 1

    fig, ax = plt.subplots(1, 1, figsize=(15, 7), dpi=300)
    # y = size in bytes
    # x = block size
    # color = codec cocktail

    colors = []

    # violin plot for kzip sizes
    violin_data_dict = defaultdict()
    for gwas_file in kzip_sizes:
        for block_size in kzip_sizes[gwas_file]:
            for codec_cocktail in kzip_sizes[gwas_file][block_size]:
                # divide size by 1e6 for readability
                try:
                    (violin_data_dict[block_size][codec_cocktail].append(
                        kzip_sizes[gwas_file][block_size][codec_cocktail] / scale))
                except KeyError:
                    try:
                        violin_data_dict[block_size][codec_cocktail] = []
                        (violin_data_dict[block_size][codec_cocktail].append(
                            kzip_sizes[gwas_file][block_size][codec_cocktail] / scale))
                    except KeyError:
                        try:
                            violin_data_dict[block_size] = {}
                            violin_data_dict[block_size][codec_cocktail] = []
                            (violin_data_dict[block_size][codec_cocktail].append(
                                kzip_sizes[gwas_file][block_size][codec_cocktail] / scale))
                        except KeyError:
                            violin_data_dict[block_size] = {}
                            violin_data_dict[block_size][codec_cocktail] = []
                            (violin_data_dict[block_size][codec_cocktail].append(
                                kzip_sizes[gwas_file][block_size][codec_cocktail] / scale))

                colors.append(codec_colors[codec_cocktail])

    violin_data = []
    for block_size in violin_data_dict:
        for codec_cocktail in violin_data_dict[block_size]:
            violin_data.append(violin_data_dict[block_size][codec_cocktail])

    x_data = []
    x_idx = 0
    for block_size in violin_data_dict:
        x_point = x_idx
        for codec_cocktail in violin_data_dict[block_size]:
            x_data.append(x_point)
            x_point += 0.2
        x_idx += 1

    # plot the violins
    plt.violinplot(violin_data,
                   x_data,
                   widths=0.15,
                   showmeans=False,
                   showmedians=True)

    # set y axis limits
    # plt.ylim(0, 4100)

    # color the violin plots by codec cocktail
    for i, pc in enumerate(plt.gca().collections):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor('black')

    # TODO: plot genomic index and pvalue index sizes
    # plot genomic index sizes
    # plot pvalue index sizes

    # FOR NOW
    # TODO: plot gzip and bgzip sizes for each file, not an average
    plt.axhline(y=gzip_sizes[0] / scale,
                color='black', linestyle='solid', linewidth=3,
                label='gzip')

    # plot the average bgzip size in blue horizontal line
    plt.axhline(y=bgzip_sizes[0] / scale,
                color='black', linestyle='dashed', linewidth=3,
                label='bgzip')



    # title
    plt.title('Compressed File Size\nby Block Size and Codec',
              fontsize=16, fontweight='bold')

    # FORMATTING
    # # log scale
    # plt.yscale('log')
    # legend, no frame
    # add custom legend for codecs
    codec_legend_elements = [Patch(facecolor=codec_colors[codec], edgecolor='black', label=codec)
                            for codec in codec_colors]

    # # add black for genomic index and grey for pvalue index
    # codec_legend_elements.append(Patch(facecolor='black', edgecolor='black', label='genomic index'))
    # codec_legend_elements.append(Patch(facecolor='grey', edgecolor='black', label='pvalue index'))

    # add bgzip dashed line to legend
    codec_legend_elements.insert(0, Line2D([0], [0],
                                           color='black', linestyle='dashed', linewidth=3,
                                           label='bgzip'))
    # add gzip solid line to legend
    codec_legend_elements.insert(1, Line2D([0], [0],
                                           color='black', linestyle='solid', linewidth=3,
                                           label='gzip'))

    # plot legend outside of plot
    plt.legend(handles=codec_legend_elements,
                title='Compression Method', frameon=False,
                loc='upper left', bbox_to_anchor=(1, 1),
                fontsize=12, title_fontsize=14)

    # label the y-axis
    plt.ylabel('Compressed file size\n(MB)', fontsize=14)

    # label the x-axis by block size
    x_ticks = range(len(violin_data_dict.keys()))
    x_tick_labels = []
    for block_size_label in violin_data_dict.keys():
        if block_size_label == 'map':
            x_tick_labels.append('map\n(1cM)')
        else:
            x_tick_labels.append(f'{block_size_label}')
    # x_tick_labels = [f'{block_size}' for block_size in violin_data_dict.keys()]
    plt.xticks(x_ticks, x_tick_labels, fontsize=12)
    plt.xlabel('Block size\n(number of records in block)', fontsize=14)

    # remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # layout
    plt.tight_layout()

    # SAVE
    plt.savefig(out_png)
