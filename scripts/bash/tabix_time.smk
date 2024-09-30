#snakemake --dryrun -s scripts/bash/tabix_time.smk --configfile data/tabix_test/tabix_config.yml
#snakemake -s scripts/bash/tabix_time.smk --configfile data/tabix_test/tabix_config.yml --rulegraph | dot -Tpdf > dag.pdf

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
        expand(f"{config.tabix_dir}{{root_file_name}}_tabix_output.txt", root_file_name=ROOT_FILE_NAMES)

rule download_bgz:
    message: "Downloading bgz files."
    input:
        manifest=f"{config.ukbb_manifest}"
    output:
        bgz_file_name=f"{config.gwas_dir}{{root_file_name}}.tsv.bgz"
    params:
        url=lambda wildcards: BGZ_URLS[ROOT_FILE_NAMES.index(wildcards.root_file_name)]
    shell:
        """
        mkdir -p {config.gwas_dir}
        cd {config.gwas_dir}
        {params.url}
        """

rule download_tabix:
    message: "Downloading tabix files."
    input:
        bgz_file_name=f"{config.gwas_dir}{{root_file_name}}.tsv.bgz"
    output:
        tbx_file_name=f"{config.gwas_dir}{{root_file_name}}.tsv.bgz.tbi"
    params:
        url=lambda wildcards: TBX_URLS[ROOT_FILE_NAMES.index(wildcards.root_file_name)]
    shell:
        """
        cd {config.gwas_dir}
        {params.url}
        """

rule tabix_search:
    message: "Searching files with tabix (applying pvalue threshold)."
    input:
        bgz_file_name=f"{config.gwas_dir}{{root_file_name}}.tsv.bgz",
        tbx_file_name=f"{config.gwas_dir}{{root_file_name}}.tsv.bgz.tbi",
        pval_indexes=f"{config.pval_indexes}"
    output:
        f"{config.tabix_dir}{{root_file_name}}_tabix_output.txt"
    shell:
        """
        mkdir -p {config.gwas_dir}
        conda activate snakemake
        cd {config.root_dir}
        python {config.scripts_dir}tabix_query.py \
        --bed {config.bed_file} \
        --gwas {input.bgz_file_name} \
        --pval_threshold {config.pval_threshold} \
        --pval_index {input.pval_indexes} \
        --out {config.tabix_dir}{wildcards.root_file_name}_tabix_output.txt
        """