import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



# how many significant SNPs are there, total?
# how many genes have significant SNPs?
def plot_snps_genes_count(trait_dict,
                          figures_dir):
    '''
    Plot the number of significant SNPs for each trait
    :param trait_dict: dictionary of traits and genes with significant SNPs
    :param figures_dir: directory to save figures
    :return: None
    '''

    # plot with two y axes, one for gene count and one for snp counts
    fig, ax1 = plt.subplots(1, 1, figsize=(10, 5))
    ax2 = ax1.twinx()

    snp_counts = {}
    gene_counts = {}
    for trait in trait_dict:
        snp_counts[trait] = 0
        trait_genes = 0

        for gene in trait_dict[trait]:
            if len(trait_dict[trait][gene]) > 0:
                trait_genes += 1
                gene_counts[trait] = trait_genes

        if trait_genes == 0:
            gene_counts[trait] = 0

        trait_genes = trait_dict[trait]
        for gene in trait_genes:
            snp_counts[trait] += len(trait_genes[gene])

    # plot gene and snp counts side by side for each trait
    x = np.arange(len(trait_dict))
    width = 0.35

    ax1.bar(x - width/2, list(gene_counts.values()), width, label='Genes',
           color='steelblue', edgecolor='black', alpha=0.7)
    ax2.bar(x + width/2, list(snp_counts.values()), width, label='SNPs',
           color='darkorange', edgecolor='black', alpha=0.7)

    # add the value of each bar at the bottom of the bar
    for i, v in enumerate(list(gene_counts.values())):
        ax1.text(i - width/2, v/2, str(v), color='black', ha='center')
    for i, v in enumerate(list(snp_counts.values())):
        ax2.text(i + width/2, v/2, str(v), color='black', ha='center')


    # labeling
    ax1.set_ylabel('Number of genes with\nsignificant SNPs', fontsize=12, color='steelblue')
    ax2.set_ylabel('Total number of\nsignificant SNPs', fontsize=12, color='darkorange')
    ax1.set_xlabel('Trait', fontsize=12)
    ax1.set_title('Gene and SNP Counts by Trait', fontsize=15, fontweight='bold')
    ax1.set_xticks(x)
    # automatically wrap text for x labels
    x_labels = list(trait_dict.keys())
    x_labels = [x_label.replace(' ', '\n') for x_label in x_labels]
    ax1.set_xticklabels(x_labels, fontsize=10)

    # one legend for both axes
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper right', frameon=False)

    # formatting
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    plt.tight_layout()
    plt.savefig(figures_dir + 'snps_genes_count.png')



# how many significant SNPs are there, per gene?
def plot_snps_per_gene(trait_dict,
                       figures_dir):
    '''
    for each train plot a separate histogram of the number of snps per gene.
    :param trait_dict: dictionary of traits and genes with significant SNPs
    :param figures_dir: directory to save figures
    :return: None
    '''

    fig, ax = plt.subplots(int(len(trait_dict)/2), 2,
                           figsize=(14, 10),
                           sharex=True, sharey=True)

    col_i = 0
    row_j = 0
    for trait in trait_dict:
        snp_counts = []
        trait_genes = trait_dict[trait]
        for gene in trait_genes:
            snp_counts.append(len(trait_genes[gene]))

        ax[row_j, col_i].hist(snp_counts, bins=20, color='darkorange', edgecolor='black', alpha=0.7)

        # log scale
        ax[row_j, col_i].set_yscale('log')

        # labeling
        ax[row_j, col_i].set_ylabel('Gene Count', fontsize=12)
        ax[row_j, col_i].set_title(trait, fontsize=15, fontweight='bold')

        # formatting
        ax[row_j, col_i].spines['top'].set_visible(False)
        ax[row_j, col_i].spines['right'].set_visible(False)

        # update subplot indices
        col_i += 1
        if col_i == 2:
            col_i = 0
            row_j += 1



    # add x-axis label to bottom subplots
    ax[-1, 0].set_xlabel('Significant SNPs per Gene', fontsize=12)
    ax[-1, 1].set_xlabel('Significant SNPs per Gene', fontsize=12)

    plt.tight_layout()
    plt.savefig(figures_dir + 'snps_per_gene.png')




# what is the distribution of p-values for significant SNPs?
def plot_pvalue_dist(trait_dict,
                     figures_dir):
    '''
    for each trait plot a histogram of the p-values of the significant SNPs
    :param trait_dict: dictionary of traits and genes with significant SNPs
    :param figures_dir: directory to save figures
    :return: None
    '''

    fig, ax = plt.subplots(int(len(trait_dict)/2), 2,
                           figsize=(14, 12),
                           sharex=True, sharey=True)

    col_i = 0
    row_j = 0
    for trait in trait_dict:
        p_values = []
        trait_genes = trait_dict[trait]
        for gene in trait_genes:
            for record in trait_genes[gene]:
                p_values.append(float(record[9]))

        ax[row_j, col_i].hist(p_values, bins=2, color='darkorange', edgecolor='black', alpha=0.7)

        # log scale
        ax[row_j, col_i].set_xscale('log')

        # labeling
        ax[row_j, col_i].set_ylabel('SNP Count', fontsize=12)
        ax[row_j, col_i].set_title(trait, fontsize=15, fontweight='bold')

        # formatting
        ax[row_j, col_i].spines['top'].set_visible(False)
        ax[row_j, col_i].spines['right'].set_visible(False)

        # update subplot indices
        col_i += 1
        if col_i == 2:
            col_i = 0
            row_j += 1

    # add x-axis label to bottom subplots
    ax[-1, 0].set_xlabel('P-Value', fontsize=12)
    ax[-1, 1].set_xlabel('P-Value', fontsize=12)

    plt.tight_layout()
    plt.savefig(figures_dir + 'pvalue_dist.png')