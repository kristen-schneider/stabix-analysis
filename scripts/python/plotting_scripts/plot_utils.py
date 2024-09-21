from collections import defaultdict

from docutils.nodes import header


def read_tabix_times(tabix_times_file):
    # gwas file: times
    tabix_times = defaultdict()
    header_1 = None
    header_2 = None
    f = open(tabix_times_file, 'r')
    for line in f:
        if len(line.strip()) == 0:
            continue
        elif "GWAS file: " in line:
            # read header1 and get gwas file name
            header_1 = line.split(',')
            gwas_file = header_1[0].strip().split(':')[1].strip()
            bytes = int(header_1[1].strip().split(':')[1].strip())
            tabix_times[gwas_file] = {}
            # read header2
            header_2 = f.readline().split(',')
        else:
            L = line.strip().split(',')
            gene = L[0]
            time = float(L[1])
            try:
                tabix_times[gwas_file][gene] = time
            except KeyError:
                tabix_times[gwas_file] = {gene: time}

    return tabix_times

def get_gene_size(bed_file):
    gene_size = dict()
    f = open(bed_file, 'r')
    for line in f:
        L = line.strip().split()
        gene = L[3]
        size = int(L[2]) - int(L[1])
        gene_size[gene] = size
    return gene_size