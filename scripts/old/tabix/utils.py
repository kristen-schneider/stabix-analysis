import tabix_utils

def read_gwas_files(gwas_list):
    # read the GWAS files from a list
    gwas_files = {}

    with open(gwas_list, 'r') as f:
        # header = f.readline()
        for line in f:
            L = line.strip().split(',')
            name = L[0]
            trait = L[1]
            samples = int(L[2])
            gwas_files[name] = (trait, samples)

    f.close()

    return gwas_files


def get_genes(bed_file):
    # {gene: {chrom: [start, end]}}
    genes = {}

    with open(bed_file, 'r') as f:
        for line in f:
            chrom, start, end, gene = line.strip().split('\t')[:4]
            if gene not in genes:
                genes[gene] = {}
            if chrom not in genes[gene]:
                genes[gene][chrom] = []
            genes[gene][chrom].append((int(start), int(end)))

    f.close()

    return genes

def check_record(record,
                 p_value_idx,
                 p_value_threshold):
    # check if the record has a p-value below the threshold
    return float(record[p_value_idx]) < p_value_threshold
