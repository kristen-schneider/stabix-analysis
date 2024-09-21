import pandas as pd
import numpy as np
import pysam
# import tabix
import time


def tabix_query_with_threshold(time_file,
                               gene_name,
                               chrm,
                               start,
                               end,
                               gwas_file,
                               p_value_idx=-1,
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

    # time tabix query
    start_time = time.time()
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

    end_time = time.time()
    elapsed_time = end_time - start_time
    time_file.write(gene_name + ',' + str(elapsed_time) + '\n')

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
    # check if the record has a p-value below the threshold
    try:
        return float(split_record[p_value_idx]) <= p_value_threshold
    except ValueError:
        return False

def read_pval_index(pval_index_file):
    pval_indexes = {}
    with open(pval_index_file, 'r') as f:
        header = f.readline()
        for line in f:
            pval = int(line.strip().split(':')[0])
            title = line.strip().split(':')[1]
            pval_indexes[pval] = title
    f.close()
    return pval_indexes










