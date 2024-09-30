from collections import defaultdict

# def read_tabix_genes(tabix_results_genes):
#     tabix_genes_times = defaultdict(dict)
#     tabix_genes_hits = defaultdict(dict)
#
#     f = open(tabix_results_genes, 'r')
#     for line in f:
#         L = line.strip().split(',')
#         if 'GWAS file: ' in L[0]:
#             gwas_file = line.strip().split(',')[0].split(':')[1].strip()
#             tabix_genes_times[gwas_file] = {}
#             tabix_genes_hits[gwas_file] = {}
#         elif 'Gene:' in L[0]:
#             gene_name, gene_time = line.split(',')
#             gene_name = gene_name.split(':')[1].strip()
#             gene_time = float(gene_time.split(':')[1].strip())
#             tabix_genes_times[gwas_file][gene_name] = gene_time
#             # if gene not in tabix_genes_hits[gwas_file], add it with value 0
#             if gene_name not in tabix_genes_hits[gwas_file]:
#                 tabix_genes_hits[gwas_file][gene_name] = 0
#         elif 'END' in L[0]:
#             gwas_file = None
#         else:
#             try:
#                 tabix_genes_hits[gwas_file][gene_name] += 1
#             except KeyError:
#                 tabix_genes_hits[gwas_file][gene_name] = 1
#
#     return tabix_genes_times, tabix_genes_hits

def read_genes(gene_results_file):
    gene_times = defaultdict()
    gene_records = defaultdict()
    with open(gene_results_file, 'r') as f:
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
                    gene_time = float(line.strip().split(',')[1].strip().split(':')[1].strip())
                    gene_times[file_name][gene] = gene_time
                    # if gene not in gene_records, add it with value 0
                    if gene not in gene_records[file_name]:
                        gene_records[file_name][gene] = 0
            else:
                try:
                    gene_records[file_name][gene] += 1
                except KeyError:
                    gene_records[file_name][gene] = 1
    f.close()

    return gene_times, gene_records