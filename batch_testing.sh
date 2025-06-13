
# setting environment
mamba create -n snakemake_env python=3.9
micromamba activate snakemake_env
micromamba install -c conda-forge snakemake statsmodels graph-tool scikit-learn multipy dill


git clone https://github.com/ayroles-lab/SBM-tools.git

cd rawData/batch
cp ../../../../Input_SBM_ab/normalized_counts_1.tsv ../../../../Input_SBM_ab/normalized_counts_2.tsv .
#
ls
# example.tsv  normalized_counts_1.tsv  normalized_counts_2.tsv




######
# collect and zip necessary files

zip example_4files.zip \
cache/blockSummary/fdr-1e-5/example/gene_block.csv \
cache/blockSummary/fdr-1e-5/example/block_summary.csv \
cache/trimmed_graph/fdr-1e-5/example.xml.gz \
cache/MCMC/blocks/fdr-1e-5/example.dill

scp chrwhe@duke.zoologi.su.se:/mnt/griffin/chrwhe/SBM_2/SBM-tools/snakemake/example_4files.zip .