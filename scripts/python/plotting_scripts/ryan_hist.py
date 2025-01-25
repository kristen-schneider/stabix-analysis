import glob
import matplotlib.pyplot as plt

def main():
    tabix_files = glob.glob('/Users/krsc0813/Library/CloudStorage/OneDrive-UCB-O365/Desktop/gwas-run/tabix_output/*')
    stabix_files = glob.glob('/Users/krsc0813/Library/CloudStorage/OneDrive-UCB-O365/Desktop/gwas-run/stabix_output/*')
    name = 'tabix'

    genes_with_sig = set()

    tabix = {}
    for tabix_file in tabix_files:
        with open(tabix_file) as f:
            line = f.readline()
            gwas = line.strip().split()[-1]
            tabix[gwas] = {}
            for line in f:
                A = line.strip().split(',')
                if len(A) == 1:
                    continue
                gene = A[0].split(':')[1][1:]
                time = float(A[1].split(':')[1][1:])
                tabix[gwas][gene] = time

    stabix = {}
    for stabix_file in stabix_files:
        with open(stabix_file) as f:
            line = f.readline()
            gwas = line.strip().split('/')[-1][:-4]
            stabix[gwas] = {}
            for line in f:
                A = line.strip().split(',')
                if len(A) == 1:
                    continue
                gene = A[0].split(':')[1][1:]
                time = int(A[1].split(':')[1][1:])/1000000
                hits = int(A[2].split(':')[1][1:])
                stabix[gwas][gene] = time
                if hits > 0:
                    genes_with_sig.add((gwas, gene))

    height = 5
    width = 8
    rows = 2
    cols = 1

    fig, axs = plt.subplots(rows, cols, figsize=(width, height), dpi=300)

    tabix_color = 'orange'
    stabix_color = 'darkblue'
    ax = axs[0]
    diffs = []
    for gwas in tabix:
        for gene in tabix[gwas]:
            if gene not in stabix[gwas]:
                print(gene)
            diffs.append( tabix[gwas][gene] / stabix[gwas][gene] )
    ax.hist(diffs, bins=100, color=stabix_color)
    # draw a dotted vertical line at x = 1
    ax.axvline(x=1, color=tabix_color, linestyle='dashed', linewidth=2)
    # align title to the left
    ax.set_title('All genes', fontsize=10, loc='left', fontweight='bold')
    ax.set_ylabel('Freq.')

    ax = axs[1]
    diffs = []
    for gwas, gene in genes_with_sig:
        if gene not in tabix[gwas]: continue
        diffs.append( tabix[gwas][gene] / stabix[gwas][gene] )
        if (tabix[gwas][gene] / stabix[gwas][gene]) < 1:
            print('tabix wins: ', gwas, gene)
    ax.hist(diffs, bins=100, color=stabix_color)
    # draw a dotted vertical line at x = 1
    ax.axvline(x=1, color=tabix_color, linestyle='dashed', linewidth=2)
    ax.set_title('Genes with significant hits', fontsize=10, loc='left', fontweight='bold')

    axs_min = min([ax.get_xlim()[0] for ax in axs])
    axs_max = max([ax.get_xlim()[1] for ax in axs])

    # super title
    fig.suptitle(f'Speedup of STABIX over {name.capitalize()}', fontsize=12, fontweight='bold')

    for ax in axs:
        #ax.set_xlim(axs_min, axs_max)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    ax.set_xlabel('Speedup')
    ax.set_ylabel('Freq.')


    plt.tight_layout()
    plt.savefig('diffs.png')

if __name__ == '__main__':
    main()