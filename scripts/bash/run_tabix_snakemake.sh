#!/usr/bin/env bash

## TODO: change SLURM options
#SBATCH -p short
#SBATCH --job-name=tabix_snakemake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=32gb
#SBATCH --time=20:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --test_output=tabix_snakemake.out
#SBATCH --error=tabix_snakemake.err

#. /opt/conda/etc/profile.d/conda.sh
#conda activate snakemake

tabix_snakemake="scripts/bash/tabix_time.smk"
config_file="test/test_tabix.yml"

#start_slice=$(date +%s.%3N)
snakemake \
    -s $tabix_snakemake \
    --cores 1 \
    --jobs 1 \
    --configfile=$config_file