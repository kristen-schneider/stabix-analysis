#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=tabix-search-gwas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=32gb
#SBATCH --time=20:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=/scratch/Users/krsc0813/bash_scripts/out/tabix-search-gwas.out
#SBATCH --error=/scratch/Users/krsc0813/bash_scripts/err/tabix-search-gwas.err

gwas_repo="/scratch/Users/krsc0813/gwas-analysis/"
gene_bed_file=$gwas_repo"data/hg19.protein_coding.bed"
gwas_files_list=$gwas_repo"data/gwas_file_names.csv"
pvalue_idx=9
pvalue_threshold=5e-8

python $gwas_repo"scripts/python/tabix/main.py"\
  --bed $gene_bed_file\
  --gwas_list $gwas_files_list\
  --gwas_dir "/scratch/Users/krsc0813/gwas_data/"\
  --figures $gwas_repo"figures/tabix/"
