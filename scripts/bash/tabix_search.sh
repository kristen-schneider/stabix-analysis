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

gene_bed_file="hg19.protein_coding.bed"
pvalue_idx=9
pvalue_threshold=5e-8

python sandbox.py
  --bed data/hg19.protein_coding.bed
  --gwas_list data/gwas_file_names.txt
  --gwas_dir data/gwas_files/
  --figures figures/tabix/