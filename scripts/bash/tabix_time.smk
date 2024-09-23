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
        f"{config.output_dir}all_gwas_data_removed.txt",
        f"{config.output_dir}all_pval_col_idx.txt",
        f"{config.output_dir}all_tabix_query_results.txt",
        f"{config.output_dir}all_tabix_query_times.txt"

rule download_bgz:
    message: "Downloading bgz files."
    input:
        manifest=f"{config.ukbb_manifest}"
    output:
        bgz_file_name=f"{config.output_dir}{{root_file_name}}.tsv.bgz"
    params:
        url=lambda wildcards: BGZ_URLS[ROOT_FILE_NAMES.index(wildcards.root_file_name)]
    shell:
        """
        mkdir -p {config.output_dir}
        cd {config.output_dir}
        {params.url}
        """

rule download_tabix:
    message: "Downloading tabix files."
    input:
        bgz_file_name=f"{config.output_dir}{{root_file_name}}.tsv.bgz"
    output:
        tbx_file_name=f"{config.output_dir}{{root_file_name}}.tsv.bgz.tbi"
    params:
        url=lambda wildcards: TBX_URLS[ROOT_FILE_NAMES.index(wildcards.root_file_name)]
    shell:
        """
        cd {config.output_dir}
        {params.url}
        """

rule get_pval_col_idx:
    message: "Getting pval column indexes."
    input:
        bgz_file_name=f"{config.output_dir}{{root_file_name}}.tsv.bgz"
    output:
        pval_file_name=f"{config.output_dir}{{root_file_name}}-pval_col_idx.txt"
    shell:
        """
        conda activate gwas_cpp
        cd {config.root_dir}
        echo GWAS file: {wildcards.root_file_name} > {output.pval_file_name}
        bgzip -dc {input.bgz_file_name} |
        head -n 1 |
        tr '\\t' '\\n' |
        grep -n -i -E 'pval|p-val|p_val' |
        cut -d ':' -f 1 -f 2 >> {output.pval_file_name}
        """

rule tabix_search:
    message: "Searching files with tabix (applying pvalue threshold)."
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
        python {config.scripts_dir}tabix_query.py \
        --bed {config.bed_file} \
        --gwas {input.bgz_file_name} \
        --pval_cols {input.pval_file_name} \
        --pval_threshold {config.pval_threshold} \
        --tabix_out {config.output_dir}{wildcards.root_file_name} \
        """

rule remove_gwas_data:
    message: "Removing GWAS files (.bgz and .tbi)."
    input:
        f"{config.output_dir}{{root_file_name}}-tabix_query_results.txt",
        f"{config.output_dir}{{root_file_name}}-tabix_query_times.txt",
    output:
        cleaned_file_name=f"{config.output_dir}{{root_file_name}}-gwas_data_removed.txt"
    params:
        bgz_file_name=lambda wildcards: f"{config.output_dir}{wildcards.root_file_name}.tsv.bgz",
        tbx_file_name=lambda wildcards: f"{config.output_dir}{wildcards.root_file_name}.tsv.bgz.tbi"
    shell:
        """
        rm {params.bgz_file_name}
        rm {params.tbx_file_name}
        touch {output.cleaned_file_name}
        """

rule remove_empty_files:
    message: "Removing empty files."
    input:
        cleaned_file_names=expand(f"{config.output_dir}{{root_file_name}}-gwas_data_removed.txt", root_file_name=ROOT_FILE_NAMES)
    output:
        f"{config.output_dir}all_gwas_data_removed.txt"
    shell:
        """
        echo {input.cleaned_file_names} | tr ' ' '\\n' > {output}
        rm {input.cleaned_file_names}
        """

rule combine_pval_col_idx:
    message: "Combining pval column indexes."
    input:
        tabix_query_results=expand(f"{config.output_dir}{{root_file_name}}-tabix_query_results.txt", root_file_name=ROOT_FILE_NAMES),
        tabix_query_times=expand(f"{config.output_dir}{{root_file_name}}-tabix_query_times.txt", root_file_name=ROOT_FILE_NAMES),
        pval_file_names=expand(f"{config.output_dir}{{root_file_name}}-pval_col_idx.txt", root_file_name=ROOT_FILE_NAMES)
    output:
        f"{config.output_dir}all_pval_col_idx.txt"
    shell:
        """
        cat {input.pval_file_names} > {output}
        rm {input.pval_file_names}
        """

rule combine_tabix_query_results:
    message: "Combining tabix query results."
    input:
        tabix_query_results=expand(f"{config.output_dir}{{root_file_name}}-tabix_query_results.txt", root_file_name=ROOT_FILE_NAMES)
    output:
        f"{config.output_dir}all_tabix_query_results.txt"
    shell:
        """
        cat {input.tabix_query_results} > {output}
        rm {input.tabix_query_results}
        """

rule combine_tabix_query_times:
    message: "Combining tabix query times."
    input:
        tabix_query_times=expand(f"{config.output_dir}{{root_file_name}}-tabix_query_times.txt", root_file_name=ROOT_FILE_NAMES)
    output:
        f"{config.output_dir}all_tabix_query_times.txt"
    shell:
        """
        cat {input.tabix_query_times} > {output}
        rm {input.tabix_query_times}
        """












