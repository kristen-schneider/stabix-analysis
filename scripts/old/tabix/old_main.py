import argparse

import plot_sig_snps as pss
import tabix_utils
import utils


num_genes = 100
p_value_threshold = 5e-8
p_value_idx = 9

def parse_args():
    parser = argparse.ArgumentParser(description='Sandbox script')
    parser.add_argument('--bed', type=str, help='gene bed file')
    parser.add_argument('--gwas_list', type=str, help='file listing GWAS files')
    parser.add_argument('--gwas_dir', type=str, help='dir with GWAS files')
    parser.add_argument('--figures', type=str, help='dir to save figures')

    return parser.parse_args()

def main():
    args = parse_args()
    bed_file = args.bed
    gwas_dir = args.gwas_dir
    gwas_files = utils.read_gwas_files(args.gwas_list)
    figures_dir = args.figures

    trait_dict = {}
    genes = utils.get_genes(bed_file)

    for gwas_file in gwas_files:

        gwas_file_path = gwas_dir + gwas_file + '.bed.gz'

        gene_count = 0
        gene_records = {}

        for gene in genes:
            gene_records[gene] = []
            for chrom in genes[gene]:
                for start, end in genes[gene][chrom]:
                    records = tabix_utils.tabix_query(chrom, start, end,
                                                      gwas_file_path)
                    gene_records[gene].extend(records)
            gene_count += 1
            if gene_count == num_genes:
                break

        trait_dict[gwas_files[gwas_file][0]] = gene_records

        x = -1

    # plotting
    pss.plot_snps_genes_count(trait_dict,
                              figures_dir)

    pss.plot_snps_per_gene(trait_dict,
                            figures_dir)

    pss.plot_pvalue_dist(trait_dict,
                            figures_dir)

    x = -1



if __name__ == '__main__':
    main()