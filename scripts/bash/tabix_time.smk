#snakemake --dryrun -s scripts/bash/tabix_time.smk --configfile data/tabix_test/tabix_config.yml -c 1 -j 1
#snakemake -s scripts/bash/tabix_time.smk --configfile data/tabix_test/tabix_config.yml --dag | dot -Tpdf > dag.pdf
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
        #expand(f"{config.output_dir}{{root_file_name}}.tsv.bgz.tbi", root_file_name=ROOT_FILE_NAMES),
        #expand(f"{config.output_dir}{{root_file_name}}-pval_col_idx.txt", root_file_name=ROOT_FILE_NAMES),
        expand(f"{config.output_dir}{{root_file_name}}-tabix_query_results.txt", root_file_name=ROOT_FILE_NAMES),
        expand(f"{config.output_dir}{{root_file_name}}-tabix_query_times.txt", root_file_name=ROOT_FILE_NAMES)

rule download_bgz:
    input:
        manifest=f"{config.ukbb_manifest}"
    output:
        bgz_file_name=f"{config.output_dir}{{root_file_name}}.tsv.bgz"
    params:
        url=lambda wildcards: BGZ_URLS[ROOT_FILE_NAMES.index(wildcards.root_file_name)]
    shell:
        """
        cd {config.output_dir}
        echo 'Downloading bgz file: {output.bgz_file_name}'
        {params.url}
        """

rule download_tabix:
    input:
        bgz_file_name=f"{config.output_dir}{{root_file_name}}.tsv.bgz"
    output:
        tbx_file_name=f"{config.output_dir}{{root_file_name}}.tsv.bgz.tbi"
    params:
        url=lambda wildcards: TBX_URLS[ROOT_FILE_NAMES.index(wildcards.root_file_name)]
    shell:
        """
        cd {config.output_dir}
        echo 'Downloading tabix file: {output.tbx_file_name}'
        {params.url}
        """

rule get_pval_col_idx:
    input:
        bgz_file_name=f"{config.output_dir}{{root_file_name}}.tsv.bgz"
    output:
        pval_file_name=f"{config.output_dir}{{root_file_name}}-pval_col_idx.txt"
    shell:
        """
        conda activate gwas_cpp
        cd {config.root_dir}
        echo 'Getting pval column index for file: {input.bgz_file_name}'
        bgzip -dc {input.bgz_file_name} |
        head -n 1 |
        tr '\\t' '\\n' |
        grep -n -i -E 'pval|p-val|p_val' |
        cut -d ':' -f 1 -f 2 > {output.pval_file_name}
        """

rule tabix_search:
    input:
        bgz_file_name=f"{config.output_dir}{{root_file_name}}.tsv.bgz",
        tbx_file_name=f"{config.output_dir}{{root_file_name}}.tsv.bgz.tbi",
        pval_file_name=f"{config.output_dir}{{root_file_name}}-pval_col_idx.txt"
    output:
        f"{config.output_dir}{{root_file_name}}-tabix_query_results.txt",
        f"{config.output_dir}{{root_file_name}}-tabix_query_times.txt"
    shell:
        """
        conda activate snakemake
        cd {config.root_dir}
        echo 'Searching file: {input.bgz_file_name}'
        python {config.scripts_dir}tabix_query.py \
        --bed {config.bed_file} \
        --gwas {input.bgz_file_name} \
        --pval_cols {input.pval_file_name} \
        --pval_threshold {config.pval_threshold} \
        --tabix_out {config.output_dir}{wildcards.root_file_name} \
        """













