shell.prefix("""
. /Users/krsc0813/miniconda3/etc/profile.d/conda.sh
conda activate snakemake;
""")

import pandas as pd
import glob
from os.path import basename

from types import SimpleNamespace
config = SimpleNamespace(**config)

# make a list of GWAS files from the manifest
PD_MANIFEST=pd.read_csv(f"{config.ukbb_manifest}")

# get list of gwas bgz files and tabix files
BGZ_FILE_NAMES = PD_MANIFEST["filename"].tolist()

# remove the .tsv.bgz extension
ROOT_FILE_NAMES = [basename(x).replace(".tsv.bgz", "") for x in BGZ_FILE_NAMES]

BGZ_URLS = PD_MANIFEST["wget"].tolist()
TBX_URLS = PD_MANIFEST["wget_tabix"].tolist()

rule all:
    input:
        expand(f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/{{root_file_name}}_2000_combo-xzb.grlz", root_file_name=ROOT_FILE_NAMES),
        expand(f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/genomic.idx", root_file_name=ROOT_FILE_NAMES),
        expand(f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/pval.idx", root_file_name=ROOT_FILE_NAMES),
        expand(f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/{{root_file_name}}.query", root_file_name=ROOT_FILE_NAMES)

rule compress_lazer:
    message: "Compressing gwas files using lazer."
    input:
        tsv_file_name=f"{config.gwas_dir}{{root_file_name}}.tsv"
    output:
        lazer_file_name=f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/{{root_file_name}}_2000_combo-xzb.grlz",
        genomic_index_file_name=f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/genomic.idx"
    shell:
        """
        ./{config.bin_dir}gwas_compress {config.config_dir}{wildcards.root_file_name}_2000_combo-xzb.yml
        """

rule index_lazer:
    message: "Indexing lazer files."
    input:
        lazer_file_name=f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/{{root_file_name}}_2000_combo-xzb.grlz",
        genomic_index_file_name=f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/genomic.idx"
    output:
        pval_index_file_name=f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/pval.idx"
    shell:
        """
        ./{config.bin_dir}gwas_index {config.config_dir}{wildcards.root_file_name}_2000_combo-xzb.yml
        """

rule decompress_lazer:
    message: "Decompressing gwas files using lazer."
    input:
        lazer_file_name=f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/{{root_file_name}}_2000_combo-xzb.grlz",
        genomic_index_file_name=f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/genomic.idx",
        pval_index_file_name=f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/pval.idx"
    output:
        tsv_file_name=f"{config.gwas_dir}{{root_file_name}}_2000_combo-xzb/{{root_file_name}}.query"
    shell:
        """
        ./{config.bin_dir}gwas_decompress {config.config_dir}{wildcards.root_file_name}_2000_combo-xzb.yml
        """