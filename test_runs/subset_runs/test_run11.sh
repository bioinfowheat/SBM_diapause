cd /mnt/griffin/chrwhe/SBM_2

run_folder=test_run11
mkdir $run_folder
cd $run_folder

git clone https://github.com/ayroles-lab/SBM-tools.git



cp /mnt/griffin/chrwhe/SBM_2/Input_SBM500/Input_SBM_ab500/normalized_counts_3.tsv /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/rawData/batch
rm /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/rawData/batch/example.tsv
cd SBM-tools/snakemake

# going to put in the parallel version of trim_networks.py to see what happens
cp /mnt/griffin/chrwhe/SBM_2/parallel_testing/annealing_SBM.py \
/mnt/griffin/chrwhe/SBM_2/parallel_testing/exportBlocks.py \
/mnt/griffin/chrwhe/SBM_2/parallel_testing/make_stacked_networks.py \
/mnt/griffin/chrwhe/SBM_2/parallel_testing/trim_networks.py \
/mnt/griffin/chrwhe/SBM_2/parallel_testing/MCMC_SBM.py /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts

# sanity check 
more /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts/make_stacked_networks.py

# setting environment
micromamba activate snakemake_env
# nano config to 
nano config.yaml
fdr: [1e-3]

# try runnning.
snakemake -s Snakefile_batchGraph --cores all &> screen.log

[Tue Jun 17 10:03:10 2025]
Error in rule exportBlocks:
    jobid: 2
    input: cache/trimmed_graph/fdr-1e-3/normalized_counts_3.xml.gz, cache/MCMC/blocks/fdr-1e-3/normalized_counts_3.dill
    output: cache/blockSummary/fdr-1e-3/normalized_counts_3/block_summary.csv, cache/blockSummary/fdr-1e-3/normalized_counts_3/gene_
block.csv
    log: logs/exportBlocks/fdr-1e-3/normalized_counts_3.log (check log file(s) for error details)

cache/
├── annealedBlock
│   └── fdr-1e-3
│       └── normalized_counts_3.dill
├── blockSummary
│   └── fdr-1e-3
│       └── normalized_counts_3
│           └── background.csv
├── equilibrate
│   └── fdr-1e-3
│       └── normalized_counts_3.dill
├── graph
│   └── normalized_counts_3.xml.gz
├── initialBlock
│   └── fdr-1e-3
│       └── normalized_counts_3.dill
├── MCMC
│   ├── blocks
│   │   └── fdr-1e-3
│   │       └── normalized_counts_3.dill
│   └── hist
│       └── fdr-1e-3
│           └── normalized_counts_3.dill
└── trimmed_graph
    └── fdr-1e-3
        └── normalized_counts_3.xml.gz


# nano config to 
nano config.yaml
fdr: [1e-2]

rm -rf cache
cd log
rm -rf *
cd ..

#### oops, deleted all.


run_folder=test_run11
mkdir $run_folder
cd $run_folder

git clone https://github.com/ayroles-lab/SBM-tools.git



cp /mnt/griffin/chrwhe/SBM_2/Input_SBM500/Input_SBM_ab500/normalized_counts_3.tsv /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/rawData/batch
rm /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/rawData/batch/example.tsv
cd SBM-tools/snakemake

# going to put in the parallel version of trim_networks.py to see what happens
cp /mnt/griffin/chrwhe/SBM_2/parallel_testing/annealing_SBM.py \
/mnt/griffin/chrwhe/SBM_2/parallel_testing/make_stacked_networks.py \
/mnt/griffin/chrwhe/SBM_2/parallel_testing/trim_networks.py \
/mnt/griffin/chrwhe/SBM_2/parallel_testing/MCMC_SBM.py /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts

# sanity check 
more /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts/make_stacked_networks.py

# setting environment
micromamba activate snakemake_env
# nano config to 
nano config.yaml
fdr: [1e-2]

# try runnning.
snakemake -s Snakefile_batchGraph --cores all &> screen.log

Error in rule exportBlocks:
    jobid: 2
    input: cache/trimmed_graph/fdr-1e-2/normalized_counts_3.xml.gz, cache/MCMC/blocks/fdr-1e-2/normalized_counts_3.dill
    output: cache/blockSummary/fdr-1e-2/normalized_counts_3/block_summary.csv, cache/blockSummary/fdr-1e-2/normalized_counts_3/gene_
block.csv
    log: logs/exportBlocks/fdr-1e-2/normalized_counts_3.log (check log file(s) for error details)

RuleException:
CalledProcessError in file /mnt/griffin/chrwhe/SBM_2/test_run11/SBM-tools/snakemake/Snakefile_batchGraph, line 131:
Command 'set -euo pipefail;  /home/chrwhe/micromamba/envs/snakemake_env/bin/python3.9 /mnt/griffin/chrwhe/SBM_2/test_run11/SBM-tools
/snakemake/.snakemake/scripts/tmpfmki3ayz.exportBlocks.py' returned non-zero exit status 1.
  File "/mnt/griffin/chrwhe/SBM_2/test_run11/SBM-tools/snakemake/Snakefile_batchGraph", line 131, in __rule_exportBlocks
  File "/home/chrwhe/micromamba/envs/snakemake_env/lib/python3.9/concurrent/futures/thread.py", line 58, in run
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: .snakemake/log/2025-06-17T104612.015006.snakemake.log
['normalized_counts_3.tsv']

# trying to fix py

cp /mnt/griffin/chrwhe/SBM_2/SBM-tools/snakemake/scripts/exportBlocks.py /mnt/griffin/chrwhe/SBM_2/test_run11/SBM-tools/snakemake/scripts

snakemake -s Snakefile_batchGraph --cores all &> screen.log

tree -L 4 cache/
cache/
├── annealedBlock
│   └── fdr-1e-2
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
│       └── normalized_counts_3.dill
├── graph
│   └── normalized_counts_3.xml.gz
├── initialBlock
│   └── fdr-1e-2
│       └── normalized_counts_3.dill
├── MCMC
│   ├── blocks
│   │   └── fdr-1e-2
│   │       └── normalized_counts_3.dill
│   └── hist
│       └── fdr-1e-2
│           └── normalized_counts_3.dill
└── trimmed_graph
    └── fdr-1e-2
        └── normalized_counts_3.xml.gz