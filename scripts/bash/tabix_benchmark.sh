#!/usr/bin/env bash

## TODO: change SLURM options
#SBATCH -p short
#SBATCH --job-name=tabix_benchmarck
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=32gb
#SBATCH --time=20:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=/scratch/Users/krsc0813/bash_scripts/out/tabix_benchmark.out
#SBATCH --error=/scratch/Users/krsc0813/bash_scripts/err/tabix_benchmark.err

#TODO: make these input options
# 0. SetUp
PANUKBB_MANIFEST="data/UKBB/PUKBB-phenotype_manifest.tsv"
PANUKBB_TRAITS_OF_INTEREST="data/UKBB/UKBB_traits.csv"
PROTEIN_CODING_BED="data/hg19.protein_coding.bed"
TEMP_DIR=""

# read the gene bed file
chrm_start_end=()
echo "Reading $PROTEIN_CODING_BED ..."
while IFS=$'\t' read -r chrm start end gene strand; do
  chrm_start_end+=($chrm:$start-$end)
done < $PROTEIN_CODING_BED


# read the traits of interest file and get a list of md5_hex for each trait
traits_of_interest_md5_hex=()
echo "Reading $PANUKBB_TRAITS_OF_INTEREST ..."
while IFS=',' read -r md5_hex_i description_i; do
  # skip the header
  if [ $md5_hex_i == "md5_hex" ]; then
    continue
  fi
  traits_of_interest_md5_hex+=($md5_hex_i)
done < $PANUKBB_TRAITS_OF_INTEREST
traits_of_interest_md5_hex+=($md5_hex_i)

#TODO: make this input options
num_GWAS=1


trait_description_col=3
trait_description_extended_col=4
number_of_cases_col=5
pops_col=11
file_name_col=26
aws_link_col=28
aws_link_tabix_col=29
wget_link_col=30
wget_link_tabix_col=31
md5_hex_col=32
size_in_byte_col=33
trait_efo_term_col=37

# start reading file
echo "Reading $PANUKBB_MANIFEST"
gwas_counter=0

while IFS= read -r trait_record; do
  # skip the header
  if [ $gwas_counter -eq 0 ]; then
    gwas_counter=$((gwas_counter+1))
    continue
  else
    # split the record by tab
    IFS=$'\t' read -r -a GWAS_fields <<< "$trait_record"

    # get the trait type
    trait_description=${GWAS_fields[$trait_description_col]}
    trait_description_extended=${GWAS_fields[$trait_description_extended_col]}
    number_of_cases=${GWAS_fields[$number_of_cases_col]}
    pops=${GWAS_fields[$pops_col]}
    file_name=${GWAS_fields[$file_name_col]}
    aws_link=${GWAS_fields[$aws_link_col]}
    aws_link_tabix=${GWAS_fields[$aws_link_tabix_col]}
    wget_link=${GWAS_fields[$wget_link_col]}
    wget_link_tabix=${GWAS_fields[$wget_link_tabix_col]}
    size_in_byte=${GWAS_fields[$size_in_byte_col]}
    md5_hex=${GWAS_fields[$md5_hex_col]}
    trait_efo_term=${GWAS_fields[$trait_efo_term_col]}

    # if looking for specific traits, search whole file for md5_hex in traits_of_interest_md5_hex
    if [ $num_GWAS -eq -1 ]; then
      if [[ " ${traits_of_interest_md5_hex[@]} " =~ " ${md5_hex} " ]]; then
        echo "GWAS File Count: $gwas_counter"
        echo "...file name:" $file_name
        echo "...trait description:" $trait_description
        echo "...trait efo term:" $trait_efo_term

        # -------------------------WORK-------------------------
        # 1. Download GWAS summary statistics
        cd $TEMP_DIR
        echo "...! downloading !..."
        echo $wget_link
#        wget $wget_link

        # 2. Search using TABIX
        echo "...! searching !..."
        # tabix search each gene
        # bed file
        echo tabix -R $PROTEIN_CODING_BED $file_name
#        tabix -R $PROTEIN_CODING_BED $file_name
#        # manual
#        for gene in ${chrm_start_end[@]}; do
#          echo "...searching for gene: $gene..."
#          # time the tabix search
#          echo time tabix $file_name $gene
##          time tabix $file_name $gene
#        done

        # -------------------------WORK-------------------------

        gwas_counter=$((gwas_counter+1))
      fi
    # if looking for the first num_GWAS traits, download the first num_GWAS GWAS summary statistics
    else
      echo "GWAS File Count: $gwas_counter"
      echo "...trait description:" $trait_description
      echo "...trait efo term:" $trait_efo_term

    # tabix search each gene
    # -------------------------WORK-------------------------
      # 1. Download GWAS summary statistics
      cd $TEMP_DIR
      echo "...! downloading !..."
      echo $wget_link
#        wget $wget_link

      # 2. Search using TABIX
      echo "...! searching !..."
      # tabix search each gene
      # bed file
      echo tabix -R $PROTEIN_CODING_BED $file_name
#        tabix -R $PROTEIN_CODING_BED $file_name
      # manual
      for gene in ${chrm_start_end[@]}; do
        echo "...searching for gene: $gene..."
        # time the tabix search
        echo time tabix $file_name $gene
#          time tabix $file_name $gene
      done
        # -------------------------WORK-------------------------

      gwas_counter=$((gwas_counter+1))
      if [ $gwas_counter -eq $((num_GWAS+1)) ]; then
        break
      fi
    fi
  fi

done < $PANUKBB_MANIFEST