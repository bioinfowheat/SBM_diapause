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


python ../../../graph_tool_sbm_hiearchy.wheat_v2.py normalized_counts_3 1e-2

scp chrwhe@duke.zoologi.su.se:/mnt/griffin/chrwhe/SBM_2/test_run9/SBM-tools/snakemake/\*.png .


rawdatafile=normalized_counts_3
pval=fdr-1e-2
zip $rawdatafile.$pval.zip \
cache/blockSummary/$pval/$rawdatafile/gene_block.csv \
cache/blockSummary/$pval/$rawdatafile/block_summary.csv \
cache/trimmed_graph/$pval/$rawdatafile.xml.gz \
cache/MCMC/blocks/$pval/$rawdatafile.dill

cp normalized_counts_3.fdr-1e-2.zip /mnt/griffin/leafar

sudo cp normalized_counts_3.fdr-1e-2.zip /mnt/│griffin/leafar                 
sudo chown leafar:leafar /mnt/griffin/leafar/normalized_counts_3.fdr-1e-2.zip                                                                                                     

 tree -L 4 cache/
cache/
├── annealedBlock
│   └── fdr-1e-2
│       ├── normalized_counts_2.dill
│       └── normalized_counts_3.dill
├── blockSummary
│   └── fdr-1e-2
│       └── normalized_counts_3
│           ├── background.csv
│           ├── block_summary.csv
│           ├── gene_block.csv
│           ├── Level_1
│           └── Level_2
├── equilibrate
│   └── fdr-1e-2
│       ├── normalized_counts_2.dill
│       └── normalized_counts_3.dill
├── graph
│   ├── normalized_counts_2.xml.gz
│   └── normalized_counts_3.xml.gz
├── initialBlock
│   └── fdr-1e-2
│       ├── normalized_counts_2.dill
│       └── normalized_counts_3.dill
├── MCMC
│   ├── blocks
│   │   └── fdr-1e-2
│   │       ├── normalized_counts_2.dill
│   │       └── normalized_counts_3.dill
│   └── hist
│       └── fdr-1e-2
│           ├── normalized_counts_2.dill
│           └── normalized_counts_3.dill
└── trimmed_graph
    └── fdr-1e-2
        ├── normalized_counts_2.xml.gz
        └── normalized_counts_3.xml.gz
        