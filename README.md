# gwas compression analysis

### Plot comparison
```
python scripts/python/plotting_scripts/plot_comparison.py
    --tabix_times data/UKBB/short_output/all_tabix_query_times.txt
    --tabix_results data/UKBB/short_output/all_tabix_query_results.txt
    --new_times data/UKBB/new_decomp_times.csv
    --new_results /not/used/yet
    --compressed_files_dir data/UKBB/new_decomp_times.csv
    --out figures/pub_figures/png/
```
<details>
  
![comparison_times](figures/pub_figures/png/tabix_vs_new_search_times.png)<br>
![comparison_sizes](figures/pub_figures/png/file_sizes.png)<br>

</details>



### Write data to be plotted
```
python scripts/python/plotting_scripts/write_tabix_data.py
   --times data/UKBB/short_output/all_tabix_query_times.txt
   --results data/UKBB/short_output/all_tabix_query_results.txt
   --bed data/bed_files/hg19.protein_coding.bed
   --out figures/pub_figures/data/short-ukbb-
```

### Plot Tabix Search Results
```
python scripts/python/plotting_scripts/plot_tabix.py
    --times figures/pub_figures/data/short-ukbb-tabix_search_times.csv
    --results figures/pub_figures/data/short-ukbb-tabix_search_results.csv
    --out figures/pub_figures/png/short-ukbb-
```

<details>
  
![tabix_times](figures/pub_figures/png/short-ukbb-tabix_search_times_hist.png)<br>
![tabix_hits](figures/pub_figures/png/short-ukbb-tabix_results_hist.png)<br>
![tabix_times_hits](figures/pub_figures/png/short-ukbb-tabix_search_times_by_gene_hits.png)<br>
![tabix_times_gene_size](figures/pub_figures/png/short-ukbb-tabix_search_times_by_gene_size_scatter.png)<br>

</details>