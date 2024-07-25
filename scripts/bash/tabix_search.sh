gene_bed_file="hg19.protein_coding.bed"
pvalue_idx=9
pvalue_threshold=5e-8

python sandbox.py
  --bed data/hg19.protein_coding.bed
  --gwas_list data/gwas_file_names.txt
  --gwas_dir data/gwas_files/
  --figures figures/tabix/