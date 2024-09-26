from collections import defaultdict
import os

def read_colors(colors_file):
    # codec, color
    colors = {}
    with open(colors_file, 'r') as f:
        for line in f:
            L = line.strip().split(',')
            colors[L[0]] = L[1]
    return colors

def read_tabix_times_data(tabix_times_file):
    # GWAS_FILE,Sig_Genes,Sig_SNPs,Tabix_Time(s)
    tabix_data = {}
    header = None
    f = open(tabix_times_file, 'r')
    for line in f:
        if header == None:
            header = line.strip().split(',')
            continue
        L = line.strip().split(',')
        gwas_file = L[0]
        genes = int(L[1])
        snps = int(L[2])
        total_tabix_time = float(L[3])

        tabix_data[gwas_file] = (genes, snps, total_tabix_time)

    return tabix_data

def read_tabix_genes(tabix_results_genes):
    tabix_genes_times = defaultdict(dict)
    tabix_genes_hits = defaultdict(dict)

    f = open(tabix_results_genes, 'r')
    for line in f:
        L = line.strip().split(',')
        if 'GWAS file: ' in L[0]:
            gwas_file = line.strip().split(',')[0].split(':')[1].strip()
            tabix_genes_times[gwas_file] = {}
            tabix_genes_hits[gwas_file] = {}
        elif 'Gene:' in L[0]:
            gene_name, gene_time = line.split(',')
            gene_name = gene_name.split(':')[1].strip()
            gene_time = float(gene_time.split(':')[1].strip())
            tabix_genes_times[gwas_file][gene_name] = gene_time
            tabix_genes_hits[gwas_file][gene_name] = 0
        else:
            tabix_genes_hits[gwas_file][gene_name] += 1
            continue

    return tabix_genes_times, tabix_genes_hits

def read_new_times_data(new_times_file):
    # GWAS_FILE,New_Time(us)
    new_data = {}
    header = None
    f = open(new_times_file, 'r')
    for line in f:
        if header == None:
            header = line.strip().split(',')
            continue
        L = line.strip().split(',')
        gwas_file = L[0]
        new_time = float(L[1])

        new_data[gwas_file] = new_time

    return new_data

def read_new_genes(new_results_genes):
    gene_times = defaultdict()
    gene_records = defaultdict()
    with open(new_results_genes, 'r') as f:
        single_file_time = 0.0
        for line in f:
            if 'GWAS' in line:
                file_name = line.strip().split(':')[1].strip().split('/')[-1].replace('.tsv', '').strip()
                gene_times[file_name] = {}
                gene_records[file_name] = {}
            elif 'Gene:' in line:
                gene = line.strip().split(',')[0].strip().split(':')[1].strip()
                # if time in line, end of gene, get time and add to single_file_time
                if 'time' in line:
                    gene_time = float(line.strip().split(',')[1].strip().split(':')[1].strip()) / 1e6
                    gene_times[file_name][gene] = gene_time
                else:
                    gene_records[file_name][gene] = 0
            else:
                gene_records[file_name][gene] += 1
                continue
    f.close()

    return gene_times, gene_records

def assign_bin(size, bins):
    for i, b in enumerate(bins):
        if size < b:
            return i - 1
    return len(bins) - 1

def get_gene_size(bed_file):
    gene_sizes = {}
    with open(bed_file, 'r') as f:
        for line in f:
            L = line.strip().split()
            gene = L[3]
            start = int(L[1])
            end = int(L[2])
            if end - start == 0:
                if gene not in gene_sizes:
                    gene_sizes[gene] = 0
                else:
                    continue
            gene_sizes[gene] = end - start
    return gene_sizes

def get_file_sizes(compressed_files_dir,
                   file_names,
                   codecs,
                   block_sizes):
    # file: (tsv, bgz, tbi, gen, pval)
    file_sizes = defaultdict(dict)
    for file in os.listdir(compressed_files_dir):
        file_name = file.replace('.tsv', '')
        if file_name in file_names:
            for block_size in block_sizes:
                for codec in codecs:
                    full_file_name = file_name + '_' + str(block_size) + '_' + codec
                    tsv_size = os.path.getsize(os.path.join(compressed_files_dir, file_name + '.tsv')) / 1e9
                    bgz_size = os.path.getsize(os.path.join(compressed_files_dir, file_name + '.tsv.bgz')) / 1e9
                    tbi_size = os.path.getsize(os.path.join(compressed_files_dir, file_name + '.tsv.bgz.tbi')) / 1e9
                    new_size = os.path.getsize(os.path.join(compressed_files_dir + full_file_name, full_file_name + '.grlz')) / 1e9
                    gen_size = os.path.getsize(os.path.join(compressed_files_dir + full_file_name, 'genomic.idx')) / 1e9
                    pval_size = os.path.getsize(os.path.join(compressed_files_dir + full_file_name, 'pval.idx')) / 1e9

                    file_sizes[full_file_name] = {'tsv': tsv_size,
                                             'bgz': bgz_size,
                                             'tbi': tbi_size,
                                             'new': new_size,
                                             'gen': gen_size,
                                             'pval': pval_size}

    return file_sizes