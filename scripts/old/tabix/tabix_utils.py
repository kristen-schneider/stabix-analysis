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
        if float(r[p_value_idx]) <= p_value_threshold:
            records.append(r)

    return records







