import argparse

import matplotlib.pyplot as plt
import os
from collections import defaultdict
import plot_utils as plot_utils

def parse_args():
    parser = argparse.ArgumentParser(description='Plot Tabix vs. XXX Results')
    parser.add_argument('--tabix_times', type=str, required=True,
                        help='tabix search times')
    parser.add_argument('--tabix_results', type=str, required=True,
                        help='tabix search results')
    parser.add_argument('--new_times', type=str, required=True,
                        help='new search times')
    parser.add_argument('--new_results', type=str, required=True,
                        help='new search results')
    parser.add_argument('--out', type=str, required=True,
                        help='output directory for plots')
    parser.add_argument('--compressed_files_dir', type=str, required=True,
                        help='directory with compressed files')
    return parser.parse_args()

def read_tabix_times_data(tabix_times_file):
    ## GWAS file: file.tsv.bgz, Bytes: 000
    ## Gene,Time(sec)
    ## xxx, 000

    # file name: sum times for all genes
    tabix_times = defaultdict(list)
    tabix_file_sizes = {}
    with open(tabix_times_file, 'r') as f:
        single_file_time = 0.0
        for line in f:
            if 'GWAS' in line:
                if single_file_time > 0:
                    try:
                        tabix_times[file_name].append(single_file_time)
                    except:
                        tabix_times[file_name] = [single_file_time]
                    single_file_time = 0.0
                file_name = line.strip().split(',')[0].strip().split(':')[1].replace('.tsv.bgz', '').strip()
                file_bytes = int(line.strip().split(',')[1].strip().split(':')[1])
                tabix_file_sizes[file_name] = file_bytes
            elif 'Gene,Time' in line:
                continue
            else:
                gene, time = line.strip().split(',')
                single_file_time += float(time)
    f.close()

    if single_file_time > 0:
        tabix_times[file_name].append(single_file_time)
    return tabix_times, tabix_file_sizes

def read_new_times_data(new_times_file):
    ## GWAS file: file.tsv.bgz, Bytes: 000
    ## task,000

    # file name: sum times for all tasks
    new_times = {}
    with open(new_times_file, 'r') as f:
        single_file_time = 0.0
        for line in f:
            if 'GWAS' in line:
                if single_file_time > 0:
                    try:
                        new_times[file_name].append(single_file_time)
                    except:
                        new_times[file_name] = [single_file_time]
                    single_file_time = 0.0
                file_name = line.strip().split(':')[1].strip().split('/')[-1].replace('.tsv', '').strip()

            else:
                task, time = line.strip().split(',')
                single_file_time += float(time.replace('μs', ''))

    f.close()

    if single_file_time > 0:
        try:
            new_times[file_name].append(single_file_time)
        except:
            new_times[file_name] = [single_file_time]


    return new_times

def get_file_size_color(file_size,
                        file_size_colors):


    color = 'NONE'
    gb_size = file_size / 1e9
    for size in file_size_colors.keys():
        if gb_size <= size:
            return file_size_colors[size]

def plot_tabix_times_hist(tabix_times,
                          tabix_file_sizes,
                          new_times,
                          out):
    # tabix times are in seconds (python time.time())
    # new times are in nanoseconds (c++ chrono::high_resolution_clock::now())


    tabix_plot_data = {}
    new_plot_data = {}

    file_size_colors = {0.5: 'purple',
                        1.0: 'blue',
                        1.5: 'green',
                        2.0: 'yellow',
                        2.5: 'orange',
                        3.0: 'red'}

    for file_name in tabix_times.keys():
        tabix_plot_data[file_name] = sum(tabix_times[file_name])
    for file_name in new_times.keys():
        new_plot_data[file_name] = sum(new_times[file_name]) / 1e9

    # x axis = tabix times
    # y axis = new times

    # plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 5), dpi=300)
    ax.set_title('Tabix vs. New Search Times', fontsize=16)
    ax.set_xlabel('Tabix Search Time (s)', fontsize=14)
    ax.set_ylabel('New Search Time (s)', fontsize=14)

    for gwas_file in new_plot_data.keys():
        try:
            color = get_file_size_color(tabix_file_sizes[gwas_file],
                                        file_size_colors)
            print(f'Plotting {gwas_file} with color {color}')
            print(f'...tabix time: {tabix_plot_data[gwas_file]}')
            print(f'...new time: {new_plot_data[gwas_file]}')
            ax.scatter(tabix_plot_data[gwas_file], new_plot_data[gwas_file],
                       marker='o',
                       s=100,
                       color=color, alpha=0.7)

        except KeyError:
            print(f'No tabix time found for {gwas_file}')

    # custom legend for file sizes
    handles = [plt.Line2D([0], [0], marker='o',
                          color='w', markerfacecolor=file_size_colors[size],
                          label=f'{size}GB', alpha=0.7)
               for size in file_size_colors.keys()]
    ax.legend(handles=handles, frameon=False, title='File Size (GB)')


    # format
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    # save
    plt.savefig(os.path.join(out, 'tabix_vs_new_search_times.png'))

def get_file_sizes(compressed_files_dir,
                   file_names):
    # file: (tsv, bgz, tbi, gen, pval)
    file_sizes = defaultdict(dict)
    for file in os.listdir(compressed_files_dir):
        file_name = file.replace('.tsv', '')
        if file_name in file_names:
            tsv_size = os.path.getsize(os.path.join(compressed_files_dir, file_name + '.tsv')) / 1e9
            bgz_size = os.path.getsize(os.path.join(compressed_files_dir, file_name + '.tsv.bgz')) / 1e9
            tbi_size = os.path.getsize(os.path.join(compressed_files_dir, file_name + '.tsv.bgz.tbi')) / 1e9
            new_size = os.path.getsize(os.path.join(compressed_files_dir + file_name + '_10000_xz/', file_name + '_10000_xz.grlz')) / 1e9
            gen_size = os.path.getsize(os.path.join(compressed_files_dir + file_name + '_10000_xz/', 'genomic.idx')) / 1e9
            pval_size = os.path.getsize(os.path.join(compressed_files_dir + file_name + '_10000_xz/', 'pval.idx')) / 1e9

            file_sizes[file_name] = {'tsv': tsv_size,
                                     'bgz': bgz_size,
                                     'tbi': tbi_size,
                                     'new': new_size,
                                     'gen': gen_size,
                                     'pval': pval_size}

    return file_sizes


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

    tabix_times, tabix_file_sizes = read_tabix_times_data(args.tabix_times)
    new_times = read_new_times_data(args.new_times)

    plot_tabix_times_hist(tabix_times,
                          tabix_file_sizes,
                          new_times,
                          args.out)

    file_names = list(tabix_times.keys())
    file_sizes = get_file_sizes(args.compressed_files_dir,
                                file_names)

    plot_file_sizes(file_sizes,
                    args.out)


if __name__ == '__main__':
    main()

