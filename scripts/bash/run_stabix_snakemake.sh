## TODO: change SLURM options
#SBATCH -p long
#SBATCH --job-name=run_xxx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=32gb
#SBATCH --time=2-23:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=out/run_xxx.out
#SBATCH --error=err/run_xxx.err

#. /opt/conda/etc/profile.d/conda.sh
#conda activate snakemake

snakemake="stabix_time-xzb.smk"
config_file="ukbb_stabix_config.yml"

snakemake\
    -s $snakemake\
    --cores 4 \
    --jobs 10 \
    --configfile=$config_file \
    --cluster-config "cluster_config.yml" \
    --cluster "sbatch -J {cluster.job-name} \\
                      -t {cluster.time} \\
                      -N {cluster.nodes} \\
                      -p {cluster.partition} \\
                      --ntasks-per-node {cluster.ntasks-per-node} \\
                      --gres={cluster.gpu} \\
                      --mem={cluster.mem} \\
                      --output {cluster.output} \\
                      --error {cluster.error}"
#    --cleanup-metadata "ten_continuous-manifest.csv"