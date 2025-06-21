cd /mnt/griffin/chrwhe/SBM_2

run_folder=test_run10
mkdir $run_folder
cd $run_folder

git clone https://github.com/ayroles-lab/SBM-tools.git



cp /mnt/griffin/chrwhe/SBM_2/Input_SBM500/Input_SBM_ab500/normalized_counts_3.tsv /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/rawData/batch
rm /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/rawData/batch/example.tsv
cd ../../snakemake

# going to put in the parallel version of trim_networks.py to see what happens
cp /mnt/griffin/chrwhe/SBM_2/parallel_testing/*.py /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts

# sanity check 
more /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts/make_stacked_networks.py

# setting environment
micromamba activate snakemake_env6
# nano config to 
nano config.yaml
fdr: [1e-2]

# try runnning.
snakemake -s Snakefile_batchGraph --cores all &> screen.log

#### crashed.