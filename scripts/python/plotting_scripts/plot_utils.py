from collections import defaultdict

from docutils.nodes import header

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

def assign_bin(size, bins):
    for i, b in enumerate(bins):
        if size < b:
            return i - 1
    return len(bins) - 1