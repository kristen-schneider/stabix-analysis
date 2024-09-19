import argparse
import tabix_utils


num_genes = 100
p_value_threshold = 5e-8
p_value_idx = 9

def parse_args():
    parser = argparse.ArgumentParser(description='Sandbox script')
    parser.add_argument('--bed', type=str, help='gene bed file')
    parser.add_argument('--gwas', type=str, help='GWAS file')

    return parser.parse_args()

def main():

    args = parse_args()

    bed_file = args.bed
    gwas_file = args.gwas

    genes = tabix_utils.get_genes(bed_file)

    gene_count = len(genes)
    

if __name__ == '__main__':
    main()
