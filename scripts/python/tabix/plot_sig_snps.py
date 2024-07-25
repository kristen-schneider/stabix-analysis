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

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    snp_counts = {}
    gene_counts = {}
    for trait in trait_dict:
        snp_counts[trait] = 0
        gene_counts[trait] = len(trait_dict[trait])
        trait_genes = trait_dict[trait]
        for gene in trait_genes:
            snp_counts[trait] += len(trait_genes[gene])

    # plot gene and snp counts side by side for each trait
    x = np.arange(len(trait_dict))
    width = 0.35

    ax.bar(x - width/2, list(gene_counts.values()), width, label='Genes')
    ax.bar(x + width/2, list(snp_counts.values()), width, label='SNPs')

    # labeling
    ax.set_ylabel('Count')
    ax.set_xlabel('Trait')
    ax.set_title('Gene and SNP Counts by Trait')
    ax.set_xticks(x)
    ax.set_xticklabels(list(trait_dict.keys()))
    ax.legend()

    # formatting
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

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

    fig, ax = plt.subplots(len(trait_dict), 1, figsize=(10, 5),
                           sharex=True, sharey=True)

    for i, trait in enumerate(trait_dict):
        snp_counts = []
        trait_genes = trait_dict[trait]
        for gene in trait_genes:
            snp_counts.append(len(trait_genes[gene]))

        ax[i].hist(snp_counts, bins=20)

        # labeling
        ax[i].set_ylabel('Gene Count')
        ax[i].set_xlabel('Number of Significant SNPs per Gene')
        ax[i].set_title(trait)

        # formatting
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)

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

    fig, ax = plt.subplots(len(trait_dict), 1, figsize=(10, 5),
                           sharex=True, sharey=True)

    for i, trait in enumerate(trait_dict):
        p_values = []
        trait_genes = trait_dict[trait]
        for gene in trait_genes:
            for record in trait_genes[gene]:
                p_values.append(float(record[9]))

        ax[i].hist(p_values, bins=20)

        # labeling
        ax[i].set_ylabel('SNP Count')
        ax[i].set_xlabel('P-Value')
        ax[i].set_title(trait)

        # formatting
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(figures_dir + 'pvalue_dist.png')