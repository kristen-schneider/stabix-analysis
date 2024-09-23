# gwas compression analysis

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
Plot tabix search times per gene, colored by bgz file size.

<details>
  
![tabix_times](figures/pub_figures/png/short-ukbb-tabix_search_times_hist.png)<br>

</details>

Plot tabix search times by gene size, colored by bgz file size.

<details>
  
![tabix_times_gene_size](figures/pub_figures/png/short-ukbb-tabix_search_times_by_gene_size_scatter.png)<br>

</details>

Plot tabix search results: number of significant hits per file.

<details>
  
![tabix_times](figures/pub_figures/png/short-ukbb-tabix_results_hist.png)<br>

</details>