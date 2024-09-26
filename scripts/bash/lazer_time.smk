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
        f"{config.output_dir}lazer_setup.done"

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

rule setup:
    message: "Setting up the environment."
    input:
        bgz_file_names=expand(f"{config.output_dir}{{root_file_name}}.tsv.bgz", root_file_name=ROOT_FILE_NAMES),
        tbx_file_names=expand(f"{config.output_dir}{{root_file_name}}.tsv.bgz.tbi", root_file_name=ROOT_FILE_NAMES)
    output:
        setup_file=f"{config.output_dir}lazer_setup.done"
    shell:
        """
        cd {config.lazer_dir}
        git submodule init
        git submodule update
        mkdir -p build
        cd build
        cmake ..
        make
        touch {output.setup_file}
        """