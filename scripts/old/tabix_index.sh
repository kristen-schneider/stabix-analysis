#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=tabix-index-gwas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=32gb
#SBATCH --time=20:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --test_output=/scratch/Users/krsc0813/bash_scripts/out/tabix-index-gwas.out
#SBATCH --error=/scratch/Users/krsc0813/bash_scripts/err/tabix-index-gwas.err

gwas_file_dir="/scratch/Users/krsc0813/gwas_data/"
gwas_file_csv="/scratch/Users/krsc0813/gwas-analysis/data/gwas_file_names.csv"

# 1. read gwas_file_csv
## header format: file_name,trait,samples,link

while IFS=, read -r gwas_file trait samples link
do
  # 1. get name of file
  echo "Processing $gwas_file"
  echo "\tTrait: $trait"
  gwas_file_path=$gwas_file_dir$gwas_file".tsv"

  ## ALL THE GWAS FILES ARE DIFF FORMAT :)
  ## 2. reformat into bed file by removing the first column and double printing the 3rd column using tab delimiter
  #echo "\t...reformatting into bed file."
  #awk -v OFS='\t' '{print $2,$3,$3,$4,$5,$6,$7,$8,$9,$10}' $gwas_file_path > $gwas_file_dir${gwas_file}.bed

  # 3. compress the file with bgzip
  echo "\t...compressing with bgzip."
  bgzip -c $gwas_file_dir${gwas_file}.bed > $gwas_file_dir${gwas_file}.bed.gz

  # 4. index the file with tabix
  echo "\t...indexing with tabix."
  tabix -p bed $gwas_file_dir${gwas_file}.bed.gz

done < $gwas_file_csv

#done < $gwas_file_list


#awk '{print $2,$3,$3,$4,$5,$6,$7,$8,$9,$10}' $gwas_file > ${gwas_file}.bed
#awk 'BEGIN{FS="\t"}{print $1,$2,$3,$5,log($11),$12}' inputfile > outputfile
#
## . bigzip and tabix index
#bgzip -c GCST90179150_buildGRCh37.bed > GCST90179150_buildGRCh37.bed.gz
#tabix -p bed GCST90179150_buildGRCh37.bed.gz
