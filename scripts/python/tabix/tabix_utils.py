import pandas as pd
import numpy as np
import pysam
# import tabix

def tabix_query_with_threshold(chrm,
                               start,
                               end,
                               gwas_file,
                               p_value_idx,
                               p_value_threshold):
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
        tb = pysam.TabixFile(gwas_file)
        # tb = tabix.open(gwas_file)
    except FileNotFoundError:
        print('File not found')
        return records

    try:
        record = tb.fetch(chrm, start, end)
        # record = tb.query(chrm, start, end)
        for r in record:
            if check_record(r, p_value_idx, p_value_threshold):
                records.append(r)
    except ValueError:
        return records

    return records

def get_genes(bed_file):
    # {gene: {chrom: [start, end]}}
    genes = {}

    with open(bed_file, 'r') as f:
        for line in f:
            chrom, start, end, gene = line.strip().split()[:4]
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
    split_record = record.split()
    # TODO: p <= 5e-08 --> p >= 7.3 (-log_10(P))
    # check if the record has a p-value below the threshold
    try:
        return float(split_record[p_value_idx]) >= p_value_threshold
    except ValueError:
        return False

def read_pval_index(pval_index_file):
    pval_indexes = {}
    with open(pval_index_file, 'r') as f:
        header = f.readline()
        for line in f:
            file_name = line.strip().split(',')[0]
            pval_index = int(line.strip().split(',')[1])
            pval_indexes[file_name] = pval_index
    f.close()
    return pval_indexes










