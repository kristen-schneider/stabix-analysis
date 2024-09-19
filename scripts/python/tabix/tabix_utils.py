import tabix

def tabix_query(chrm,
                start,
                end,
                gwas_file,
                p_value_idx=9,
                p_value_threshold=5e-8):
    '''
    Query a GWAS file for one gene using tabix
    :param chrm: chromosome
    :param start: bp start
    :param end: bp end
    :param gwas_file: GWAS file
    :param p_value_idx: index of the p-value column
    :param p_value_threshold: p-value threshold
    :return: list of records at or below the p-value threshold for the gene
    '''
    records = []

    try:
        tb = tabix.open(gwas_file)
    except tabix.TabixError:
        print('Error opening tabix file')
        return records

    record = tb.query(chrm, start, end)
    for r in record:
        if check_record(r, p_value_idx, p_value_threshold):
            records.append(r)

    return records

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







