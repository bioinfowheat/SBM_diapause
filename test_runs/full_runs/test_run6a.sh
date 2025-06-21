# restarting run6, as it seems to be hung on 
tail *log
niter:   935  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1381973.6  ΔS:      33.5061  moves:  5134 
niter:   936  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1381872.6  ΔS:     -100.953  moves:  4247 
niter:   937  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1382043.9  ΔS:      171.329  moves:  4799 
niter:   938  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1381941.2  ΔS:     -102.741  moves:  5056 


# snakemake_env2 and snakemake_env3 are free.

cd /mnt/griffin/chrwhe/SBM_2
cp -r test_run6 test_run6a

cd /mnt/griffin/chrwhe/SBM_2/test_run6a/SBM-tools/snakemake
rm screen.log
micromamba activate snakemake_env2


tree -L 4 cache
cache
├── annealedBlock
│   └── fdr-1e-2
│       └── normalized_counts.dill
├── equilibrate
│   └── fdr-1e-2
│       └── normalized_counts.dill
├── graph
│   └── normalized_counts.xml.gz
├── initialBlock
│   └── fdr-1e-2
│       └── normalized_counts.dill
├── MCMC
│   ├── blocks
│   │   └── fdr-1e-2
│   └── hist
│       └── fdr-1e-2
└── trimmed_graph
    └── fdr-1e-2
        └── normalized_counts.xml.gz

# restarting
snakemake -s Snakefile_batchGraph --cores all --unlock argument &> screen.log
snakemake -s Snakefile_batchGraph --cores all &> screen.log

# huh ... looks like it needs to go back and start at trim ... 
Job stats:
job             count
------------  -------
GO                  1
MCMC                1
all                 1
annealing           1
equilibrate         1
exportBlocks        1
minSBM              1
trim                1
total               8

Select jobs to execute...

[Mon Jun 16 14:05:28 2025]
rule trim:
    input: cache/graph/normalized_counts.xml.gz
    output: cache/trimmed_graph/fdr-1e-2/normalized_counts.xml.gz
    log: logs/trim/fdr-1e-2/normalized_counts.log
    jobid: 3
    reason: Updated input files: cache/graph/normalized_counts.xml.gz
    wildcards: fdr=1e-2, label=normalized_counts
    resources: tmpdir=/tmp

# stuck here for days
 tail *log
niter:   929  count:    0  breaks:  0  min_S: 1380983.6  max_S: 1381995.0  S: 1381235.8  ΔS:      30.8847  moves:  6465 
niter:   930  count:    0  breaks:  0  min_S: 1380983.6  max_S: 1381995.0  S: 1381302.2  ΔS:      66.4583  moves:  4993 
niter:   931  count:    0  breaks:  0  min_S: 1380983.6  max_S: 1381995.0  S: 1381354.2  ΔS:      51.9377  moves:  3969 
niter:   932  count:    0  breaks:  0  min_S: 1380983.6  max_S: 1381995.0  S: 1381298.8  ΔS:     -55.3314  moves:  4571 
niter:   933  count:    0  breaks:  0  min_S: 1380983.6  max_S: 1381995.0  S: 1381199.3  ΔS:     -99.5010  moves:  3979 
niter:   934  count:    0  breaks:  0  min_S: 1380983.6  max_S: 1381995.0  S: 1381220.1  ΔS:      20.8134  moves:  3328 
niter:   935  count:    0  breaks:  0  min_S: 1380983.6  max_S: 1381995.0  S: 1381250.3  ΔS:      30.1617  moves:  5432 
niter:   936  count:    0  breaks:  0  min_S: 1380983.6  max_S: 1381995.0  S: 1381127.8  ΔS:     -122.540  moves:  5517 
niter:   937  count:    0  breaks:  0  min_S: 1380983.6  max_S: 1381995.0  S: 1381069.1  ΔS:     -58.6153  moves:  6091 
niter:   938  count:    0  breaks:  0  min_S: 1380983.6  max_S: 1381995.0  S: 1381145.6  ΔS:      76.4102  moves:  5009 

######
######
# restarting

run_folder=test_run6a
cd /mnt/griffin/chrwhe/SBM_2/$run_folder
rm -rf SBM*

git clone https://github.com/ayroles-lab/SBM-tools.git



cd SBM-tools/rawData/batch
cp /mnt/griffin/chrwhe/SBM_2/Input_SBM_ab/normalized_counts.tsv .
rm example*  
cd ../../snakemake

# going to put in the parallel version of trim_networks.py to see what happens
cp /mnt/griffin/chrwhe/SBM_2/parallel_testing/annealing_SBM.py \
/mnt/griffin/chrwhe/SBM_2/parallel_testing/make_stacked_networks.py \
/mnt/griffin/chrwhe/SBM_2/parallel_testing/trim_networks.py \
/mnt/griffin/chrwhe/SBM_2/parallel_testing/MCMC_SBM.py /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts


# sanity check 
more /mnt/griffin/chrwhe/SBM_2/test_run6a/SBM-tools/snakemake/scripts/make_stacked_networks.py

# setting environment
micromamba activate snakemake_env6
# nano config to 
nano config.yaml
fdr: [1e-2]

# try runnning.
snakemake -s Snakefile_batchGraph --cores all &> screen.log


tree -L 4 cache/                                  
├── annealedBlock
│   └── fdr-1e-2
│       └── normalized_counts.dill
├── blockSummary
│   └── fdr-1e-2
│       └── normalized_counts
│           ├── background.csv
│           ├── block_summary.csv
│           ├── gene_block.csv
│           ├── Level_1
│           ├── Level_2
│           ├── Level_4
│           └── Level_5
├── equilibrate
│   └── fdr-1e-2
│       └── normalized_counts.dill
├── graph
│   └── normalized_counts.xml.gz
├── initialBlock
│   └── fdr-1e-2
│       └── normalized_counts.dill
├── MCMC
│   ├── blocks
│   │   └── fdr-1e-2
│   │       └── normalized_counts.dill
│   └── hist
│       └── fdr-1e-2
│           └── normalized_counts.dill
└── trimmed_graph
    └── fdr-1e-2
        └── normalized_counts.xml.gz


###########
# plotting
cd /mnt/griffin/chrwhe/SBM_2/test_run6a/SBM-tools/snakemake
# getting output files for local plotting.
rawdatafile=normalized_counts
pval=fdr-1e-2
zip $rawdatafile.$pval.zip \
cache/blockSummary/$pval/$rawdatafile/gene_block.csv \
cache/blockSummary/$pval/$rawdatafile/block_summary.csv \
cache/trimmed_graph/$pval/$rawdatafile.xml.gz \
cache/MCMC/blocks/$pval/$rawdatafile.dill


sudo cp normalized_counts.fdr-1e-2.zip /mnt/│griffin/leafar                 
sudo chown leafar:leafar /mnt/griffin/leafar/normalized_counts.fdr-1e-2.zip                                                                                                    

