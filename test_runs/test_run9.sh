cd /mnt/griffin/chrwhe/SBM_2

run_folder=test_run9
mkdir $run_folder
cd $run_folder

git clone https://github.com/ayroles-lab/SBM-tools.git


cd $run_folder/SBM-tools/rawData/batch
rm example*  

cp /mnt/griffin/chrwhe/SBM_2/Input_SBM500/Input_SBM_ab500/normalized_counts_2.tsv /mnt/griffin/chrwhe/SBM_2/Input_SBM500/Input_SBM_ab500/normalized_counts_3.tsv /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/rawData/batch
cd /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/rawData/batch
cd ../../snakemake

# going to put in the parallel version of trim_networks.py to see what happens
cp /mnt/griffin/chrwhe/SBM_2/parallel_testing/trim_networks.py /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts
cp /mnt/griffin/chrwhe/SBM_2/parallel_testing/MCMC_SBM.py /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts
cp /mnt/griffin/chrwhe/SBM_2/parallel_testing/make_stacked_networks.py /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts

# sanity check
cd /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts
more make_stacked_networks.py

# setting environment
micromamba activate snakemake_env6
# nano config to 
nano config.yaml
fdr: [1e-2]

# try runnning.
snakemake -s Snakefile_batchGraph --cores all &> screen.log
