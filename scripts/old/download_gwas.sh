#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=download-gwas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=32gb
#SBATCH --time=20:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --test_output=/scratch/Users/krsc0813/bash_scripts/out/download-gwas.out
#SBATCH --error=/scratch/Users/krsc0813/bash_scripts/err/download-gwas.err

# download many GWAS summary statistics
# http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/

data_dir="/scratch/Users/krsc0813/gwas_data/"
cd $data_dir

gwas_url_list=(
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179114/GCST90179114_buildGRCh37.tsv"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179115/GCST90179115_buildGRCh37.tsv"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179117/GCST90179117_buildGRCh37.tsv"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179149/GCST90179149_buildGRCh37.tsv"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179150/GCST90179150_buildGRCh37.tsv"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179151/GCST90179151_buildGRCh37.tsv"
    )

for url in "${gwas_url_list[@]}"
do
    echo $url
    wget $url
done
