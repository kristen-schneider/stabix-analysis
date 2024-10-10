# gwas compression analysis

### Run tabix experiments (snakemake)
[stabix_mamba.yml](https://github.com/kristen-schneider/stabix/blob/main/stabix_mamba.yml)
```
mamba create -n stabix -f stabix_mamba.yml
mamba activate stabix
snakemake -s scripts/bash/tabix_time.smk
```
or
```
mamba create -n stabix -f stabix_mamba.yml
mamba activate stabix
sbatch run_tabix_snakemake.sh
```

### Run STABIX experiments (snakemake)
```
mamba create -n stabix -f stabix_mamba.yml
mamba activate stabix
snakemake -s scripts/bash/stabix_time.smk
```
or
```
mamba create -n stabix -f stabix_mamba.yml
mamba activate stabix
sbatch run_stabix_snakemake.sh
```

### Plot main hexagonal heatmap
```
python scripts/python/plotting_scripts/publish.py \
 --data /path/to/data/ \
 --bed data/bed_files/hg19.protein_coding.bed \
 --out figures/
```

<details>

![hex_heatmap](figures/gene_times.png)<br>
![file_sizes](figures/file_sizes_bar.png)<br>
![sepeedup](scripts/python/plotting_scripts/diffs.png)<br>

</details>

| file sizes                              | num blocks                              |
|-----------------------------------------|-----------------------------------------|
| ![stabix_sizes](figures/file_sizes.png) | ![stabix_sizes](figures/num_blocks.png) |
| p-value hits                            | sig SNPs                                |
| ![stabix_sizes](figures/pval_hits.png)  | ![stabix_sizes](figures/snps_hits.png)  |


### Plot combination results
```
python scripts/python/plotting_scripts/plot_compare.py \
 --root /path/to/data/ \
 --bed data/bed_files/hg19.protein_coding.bed \
 --out figures/
```

<details>

![stabix_single](figures/continuous-103220-both_sexes_compare_times.png)<br>

</details>

### Plot column results
```
python scripts/python/plotting_scripts/plot_columns.py \
 --decomp /path/to/data/ \
 --colors figures/colors.txt \
 --out figures/
   
```
<details>
  
![columns](figures/column_decompression.png)<br>

</details>


