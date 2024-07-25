#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=download-gwas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=32gb
#SBATCH --time=20:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=/scratch/Users/krsc0813/bash_scripts/out/download-gwas.out
#SBATCH --error=/scratch/Users/krsc0813/bash_scripts/err/download-gwas.err

# download many GWAS summary statistics
# http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/

data_dir="/scratch/Users/krsc0813/gwas_data/"
cd $data_dir

gwas_url_list=(
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179150/GCST90179150_buildGRCh37.tsv"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179151/GCST90179151_buildGRCh37.tsv"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179151/GCST90179152_buildGRCh37.tsv"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179179/GCST90179179_buildGRCh37.tsv.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179320/GCST90179320_buildGRCh38.tsv.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90179001-GCST90180000/GCST90179411/GCST90179411_buildGRCh38.tsv.gz"
    )

for url in "${gwas_url_list[@]}"
do
    echo $url
    wget $url
done
