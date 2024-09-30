from collections import defaultdict
import os

from docutils.nodes import header


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

def read_genes(gene_results_file,
               micro):
    gene_times = defaultdict()
    gene_records = defaultdict()
    gene_pval_hits = defaultdict()
    with open(gene_results_file, 'r') as f:
        single_file_time = 0.0
        for line in f:
            if 'GWAS' in line:
                file_name = line.strip().split(':')[1].strip().split('/')[-1].replace('.tsv', '').strip()
                gene_times[file_name] = {}
                gene_records[file_name] = {}
                gene_pval_hits[file_name] = {}
            elif 'Gene:' in line:
                gene = line.strip().split(',')[0].strip().split(':')[1].strip()
                # if time in line, end of gene, get time and add to single_file_time
                if 'time' in line:
                    if gene == 'SLC39A9':
                        x = 1
                    if micro:
                        gene_time = float(line.strip().split(',')[1].strip().split(':')[1].strip()) / 1e6
                        pval_hit = int(line.strip().split(',')[2].strip().split(':')[1].strip())
                    else:
                        gene_time = float(line.strip().split(',')[1].strip().split(':')[1].strip())
                        pval_hit = -1
                    try:
                        gene_times[file_name][gene] += gene_time
                        gene_pval_hits[file_name][gene] += pval_hit
                    except KeyError:
                        gene_times[file_name][gene] = gene_time
                        gene_pval_hits[file_name][gene] = pval_hit
                    # if gene not in gene_records, add it with value 0
                    if gene not in gene_records[file_name]:
                        gene_records[file_name][gene] = 0
            else:
                try:
                    gene_records[file_name][gene] += 1
                except KeyError:
                    gene_records[file_name][gene] = 1
    f.close()
    return gene_times, gene_records, gene_pval_hits

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
            try:
                gene_sizes[gene] += end - start
            except KeyError:
                gene_sizes[gene] = end - start
    return gene_sizes

# get gene_instances from bed file
def get_gene_instances(bed_file):
    gene_instances = {}
    with open(bed_file, 'r') as f:
        for line in f:
            L = line.strip().split()
            gene = L[3]
            try:
                gene_instances[gene] += 1
            except KeyError:
                gene_instances[gene] = 1
    return gene_instances

def get_file_sizes(compressed_files_dir,
                   file_names,
                   codecs,
                   block_sizes):
    # file: (tsv, bgz, tbi, gen, pval)
    file_sizes = defaultdict(dict)
    for file in os.listdir(compressed_files_dir):
        file_name = file.replace('.tsv', '')
        if file_name in file_names:
            file_sizes[file_name] = {}
            for block_size in block_sizes:
                for codec in codecs:
                    full_file_name = file_name + '_' + str(block_size) + '_' + codec
                    tsv_size = os.path.getsize(os.path.join(compressed_files_dir, file_name + '.tsv')) / 1e9
                    bgz_size = os.path.getsize(os.path.join(compressed_files_dir, file_name + '.tsv.bgz')) / 1e9
                    tbi_size = os.path.getsize(os.path.join(compressed_files_dir, file_name + '.tsv.bgz.tbi')) / 1e9
                    new_size = os.path.getsize(os.path.join(compressed_files_dir + full_file_name, full_file_name + '.grlz')) / 1e9
                    gen_size = os.path.getsize(os.path.join(compressed_files_dir + full_file_name, 'genomic.idx')) / 1e9
                    pval_size = os.path.getsize(os.path.join(compressed_files_dir + full_file_name, 'pval.idx')) / 1e9

                    try:
                        file_sizes[file_name][block_size][codec] = {'tsv': tsv_size,
                                                                    'bgz': bgz_size,
                                                                    'tbi': tbi_size,
                                                                    'new': new_size,
                                                                    'gen': gen_size,
                                                                    'pval': pval_size}
                    except KeyError:
                        file_sizes[file_name][block_size] = {codec: {'tsv': tsv_size,
                                                                    'bgz': bgz_size,
                                                                    'tbi': tbi_size,
                                                                    'new': new_size,
                                                                    'gen': gen_size,
                                                                    'pval': pval_size}}
    return file_sizes

def check_gene_records(tabix_gene_records,
                       lzr_gene_records):

    # assert that there are the same number of genes in both records
    assert len(tabix_gene_records) == len(lzr_gene_records)

    for gene in tabix_gene_records.keys():
        if lzr_gene_records[gene] != tabix_gene_records[gene]:
            print('Gene:', gene)
            print('--LZR:', lzr_gene_records[gene])
            print('--Tabix:', tabix_gene_records[gene])

def read_bed_file(bed_file):
    bed_dict = {}
    with open(bed_file, 'r') as f:
        for line in f:
            L = line.strip().split()
            try:
                chrm = int(L[0])
            except ValueError:
                if L[0] == 'X':
                    chrm = 23
                elif L[0] == 'Y':
                    chrm = 24
                elif L[0] == 'MT':
                    chrm = 25
            bp_start = int(L[1])
            bp_end = int(L[2])
            gene = L[3]
            try:
                bed_dict[gene].append((chrm, bp_start, bp_end))
            except KeyError:
                bed_dict[gene] = [(chrm, bp_start, bp_end)]
    return bed_dict

def read_genomic_index(genomic_index_file):
    chrm_index = {}
    with open(genomic_index_file, 'r') as f:
        header = f.readline()
        for line in f:
            L = line.strip().split(',')
            block_ID = int(L[0])
            chrm = int(L[1])
            bp_start = int(L[2])
            try:
                chrm_index[chrm].append(bp_start)
            except KeyError:
                chrm_index[chrm] = [bp_start]

    return chrm_index

def get_bp_index_from_genomic_index(genomic_index,
                                   gene_locations):
    # get the index of the base pair in the genomic index
    # get the list of base pairs for the chromosome
    gene_blocks = {}
    for gene in gene_locations:
        # if gene == 'NLGN4Y':
        #     x = gene_locations[gene]
        #     z = 1
        single_gene_blocks = 0
        for location in gene_locations[gene]:
            chrm = location[0]
            bp_start = location[1]
            bp_end = location[2]
            start_idx = None
            end_idx = None

            # if the chromosome is Y, or MT (for current file these are not in the gwas)
            if chrm == 24 or chrm == 25:
                break
            # find the index of the start base pair in the genomic index
            # then find the index of the end base pair in the genomic index
            # then add the indexes to the list

            for i, bp, in enumerate(genomic_index[chrm]):
                if bp > bp_start:
                    start_idx = i - 1
                    break
            if start_idx == None:
                # if chrm 1, start index = 0
                if chrm == 1:
                    start_idx = 0
                    break
                # if chrm > 1, start index = previous chrm end index
                else:
                    start_idx = -1
                    break

            for i, bp in enumerate(genomic_index[chrm]):
                if bp > bp_end:
                    end_idx = i - 1
                    break
            if end_idx == None:
                end_idx = len(genomic_index[chrm]) - 1

            single_gene_blocks += end_idx - start_idx + 1
        gene_blocks[gene] = single_gene_blocks
    return gene_blocks

def read_column_decompression(decompression_times_file,
                              block_size,
                              decom_times,
                              comp_sizes):

    # col_idx_data_type = {0: 'int', 1: 'int',
    #                      2: 'string', 3: 'string', 8: 'string',
    #                      4: 'float', 5: 'float', 6: 'float', 7: 'float'}

    # 0 and 1 = int
    ints = [0, 1]
    # 2, 3, 38-43 = string
    strings = [2, 3]
    strings += list(range(38, 44))
    # 4-37 = float
    floats = list(range(4, 38))
    col_idx_data_type = {}
    for i in ints:
        col_idx_data_type[i] = 'int'
    for i in strings:
        col_idx_data_type[i] = 'string'
    for i in floats:
        col_idx_data_type[i] = 'float'




    with open(decompression_times_file, 'r') as f:
        # read 2 headers
        f.readline()
        f.readline()
        for line in f:
            L = line.strip().split(',')
            column_idx = L[1]
            time = int(L[2])
            size = int(L[3])
            codec = L[4]
            data_type = col_idx_data_type[int(column_idx)]

            try:
                decom_times[block_size][codec][data_type].append(time)
                comp_sizes[block_size][codec][data_type].append(size)
            except KeyError:
                try:
                    decom_times[block_size][codec][data_type] = [time]
                    comp_sizes[block_size][codec][data_type] = [size]
                except KeyError:
                    try:
                        decom_times[block_size][codec] = {data_type: [time]}
                        comp_sizes[block_size][codec] = {data_type: [size]}
                    except KeyError:
                        decom_times[block_size] = {codec: {data_type: [time]}}
                        comp_sizes[block_size] = {codec: {data_type: [size]}}

    return decom_times, comp_sizes
