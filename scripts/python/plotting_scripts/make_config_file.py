import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(description='Write table of info for files')
    parser.add_argument('--gwas', type=str, required=True,
                        help='dir with data')
    parser.add_argument('--pvals', type=str, required=True,
                        help='pval column idx for each gwas')
    parser.add_argument('--out', type=str, required=True,
                        help='config directory for output data')
    return parser.parse_args()


def read_pvals(pvals_file):
    pvals = {}
    with open(pvals_file, 'r') as f:
        header = f.readline()
        for line in f:
            gwas, pval = line.strip().split(',')
            pvals[gwas] = float(pval)
    return pvals


def write_config_file(gwas,
                      defaults,
                      gwas_pval,
                      out):
    gwas_basename = os.path.basename(gwas).replace('.tsv', '')
    config_file = out + gwas_basename + '.yml'
    with open(config_file, 'w') as f:
        f.write('gwas_file:\n\t{}\n\n'.format(gwas))
        f.write('block_size:\n\t{}\n\n'.format(defaults['block_size']))
        f.write('out_name:\n\t{}\n\n'.format(defaults['out_name']))
        f.write('query:\n')
        f.write('\tgenomic:\n\t\t{}\n'.format(defaults['query_genomic']))
        f.write('\tcol_idx:\n\t\t{}\n'.format(gwas_pval))
        f.write('\tbins:\n\t\t{}\n'.format(defaults['bins']))
        f.write('\tthreshold:\n\t\t{}\n\n'.format(defaults['threshold']))
        f.write('extra_index:\n\t{}\n\n'.format(defaults['extra_index']))
        f.write('codecs:\n')
        f.write('\tint:\n\t\t{}\n'.format(defaults['codecs_int']))
        f.write('\tfloat:\n\t\t{}\n'.format(defaults['codecs_float']))
        f.write('\tstring:\n\t\t{}\n'.format(defaults['codecs_string']))

    f.close()

def main():
    defaults = {
        'out_name': 'combo-xzb',
        'block_size': '2000',
        'query_genomic': 'hg19.protein_coding.bed',
        'bins': '0.3, 4.3, 7.29',
        'threshold': '>= 7.3',
        'extra_index': 'pval',
        'codecs_int': 'xz',
        'codecs_float': 'zlib',
        'codecs_string': 'bz2',
    }

    args = parse_args()
    gwas = args.gwas
    gwas_basename = os.path.basename(gwas).replace('.tsv', '')
    pvals = read_pvals(args.pvals)
    gwas_pval = pvals[gwas_basename]
    out = args.out

    write_config_file(gwas,
                        defaults,
                        gwas_pval,
                        out)

if __name__ == '__main__':
    main()
