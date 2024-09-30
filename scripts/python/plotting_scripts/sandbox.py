import matplotlib.pyplot as plt
import numpy as np

import plot_utils as pltu

bed_file = '/Users/krsc0813/PycharmProjects/gwas-analysis/data/bed_files/hg19.protein_coding.bed'
idx_2000 = '/Users/krsc0813/PycharmProjects/gwas-analysis/data/UKBB/continuous-103220-both_sexes_output/zlib_2000.gen.idx'

genes = pltu.read_bed_file(bed_file)
idx = pltu.read_genomic_index(idx_2000)

fig, ax = plt.subplots(1, 1, figsize=(30, 10), dpi=300)

for chrm in idx.keys():
    # y = chrm
    # x = start bp to end bp
    y = chrm
    x_0 = 0
    x_1 = idx[chrm][-1]
    ax.plot([x_0, x_1], [y, y],
            color='black',
            linewidth=2.5)

num_genes = 0
for gene in genes:
    for loc in genes[gene]:
        y = loc[0]
        x_0 = loc[1]
        x_1 = loc[2]
        ax.plot([x_0, x_1], [y, y],
                color='red',
                linewidth=2.5)
    num_genes += 1
    if num_genes > 8000:
        break

ax.set_yticks(list(idx.keys()))
ax.set_yticklabels(list(idx.keys()))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

plt.tight_layout()
plt.savefig('/Users/krsc0813/PycharmProjects/gwas-analysis/figures/pub_figures/sandbox.png')


# for each chromosome in bed file plot a horizontal line for start bp to end bp
# for each gene in the bed file, plot a horizonal line for start bp to end bp
