from typing import List


def get_column_names(tsv_file: str) -> List[str]:
    with open(tsv_file, 'r') as f:
        header = next(f).strip().split()
        # 'chr' -> 'chrm'
        header = list(map(lambda x: 'chrm' if x == 'chr' else x, header))
    return header


def get_genes(bed_file):
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
