gwas_file_dir="data/gwas_files/"
gwas_file_csv="data/gwas_files.csv"

# 1. read gwas_file_csv
## header format: file_name,trait,samples,link

while IFS=, read -r gwas_file trait samples link
do
  # 1. get name of file
  echo "Processing $gwas_file"
  echo "\tTrait: $trait"
  gwas_file_path=$gwas_file_dir$gwas_file".tsv"

  # 2. reformat into bed file by removing the first column and double printing the 3rd column using tab delimiter
  echo "\t...reformatting into bed file."
  awk -v OFS='\t' '{print $2,$3,$3,$4,$5,$6,$7,$8,$9,$10}' $gwas_file_path > $gwas_file_dir${gwas_file}.bed

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