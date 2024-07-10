# gwas compression analysis

### Plot file sizes by codec and block size
```
python scripts/plot_file_sizes.py \
    --file_names data/gwas_file_names.txt \
    --file_sizes data/gwas_file_sizes.csv\ 
    --colors figures/colors.txt \
    --plot_dir figures/
```
<details>
  
![file_sizes](figures/GCST90179150_buildGRCh37_file_sizes.png)<br>

</details>

### Plot column sizes by codec and block size
```
python src/plot_col_sizes.py \
    --col_size_dir data/col_sizes \
     --output_dir figures/ \
     --colors figures/colors.txt
```
<details>

| block size | with fastpfor                             | without fastpfor                            |
|------------|-------------------------------------------|---------------------------------------------|
| 2000       | ![2000](figures/col_sizes_2000_fpf.png)   | ![2000](figures/col_sizes_2000_nofpf.png)   |
| 5000       | ![5000](figures/col_sizes_5000_fpf.png)   | ![5000](figures/col_sizes_5000_nofpf.png)   |
| 10000      | ![10000](figures/col_sizes_10000_fpf.png) | ![10000](figures/col_sizes_10000_nofpf.png) |
| 20000      | ![20000](figures/col_sizes_20000_fpf.png) | ![20000](figures/col_sizes_20000_nofpf.png) |

</details>