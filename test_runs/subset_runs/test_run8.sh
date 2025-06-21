# trying to start a new run with some of the generated files from a previous 
# run. copying those into folders to see if it detects them and starts where they left off

# grap new test dataset
scp Input_SBM500.zip chrwhe@duke.zoologi.su.se:/mnt/griffin/chrwhe/SBM_2/
unzip Input_SBM500.zip




cd /mnt/griffin/chrwhe/SBM_2

micromamba deactivate


cd /mnt/griffin/chrwhe/SBM_2
mkdir test_run8
cd test_run8

git clone https://github.com/ayroles-lab/SBM-tools.git


cd SBM-tools/rawData/batch
rm example*  

cp /mnt/griffin/chrwhe/SBM_2/Input_SBM500/Input_SBM_ab500/normalized_counts_1.tsv /mnt/griffin/chrwhe/SBM_2/test_run8/SBM-tools/rawData/batch
cd /mnt/griffin/chrwhe/SBM_2/test_run8/SBM-tools/rawData/batch
cd ../../snakemake

# going to put in the parallel version of trim_networks.py to see what happens
cp /mnt/griffin/chrwhe/SBM_2/parallel_testing/trim_networks.py /mnt/griffin/chrwhe/SBM_2/test_run8/SBM-tools/snakemake/scripts
# sanity check
cd /mnt/griffin/chrwhe/SBM_2/test_run8/SBM-tools/snakemake/scripts
more trim_networks.py

# setting environment
micromamba activate snakemake_env6
# nano config to 
fdr: [1e-2]

# try runnning.
snakemake -s Snakefile_batchGraph --cores all &> screen.log


[Mon Jun 16 17:13:00 2025]
Finished job 6.
5 of 9 steps (56%) done
Select jobs to execute...

[Mon Jun 16 17:13:00 2025]
rule MCMC:
    input: cache/trimmed_graph/fdr-1e-2/normalized_counts_1.xml.gz, cache/equilibrate/fdr-1e-2/normalized_counts_1.dill
    output: cache/MCMC/blocks/fdr-1e-2/normalized_counts_1.dill, cache/MCMC/hist/fdr-1e-2/normalized_counts_1.dill
    log: logs/MCMC/fdr-1e-2/normalized_counts_1.log
    jobid: 5
    reason: Missing output files: cache/MCMC/blocks/fdr-1e-2/normalized_counts_1.dill; Input files updated by another job: cache/tri
mmed_graph/fdr-1e-2/normalized_counts_1.xml.gz, cache/equilibrate/fdr-1e-2/normalized_counts_1.dill
    wildcards: fdr=1e-2, label=normalized_counts_1
    threads: 16
    resources: tmpdir=/tmp

niter:     1  count:    0  breaks:  0  min_S: 6889.3473  max_S: 6891.3795  S: 6889.3473  ΔS:     -2.03227  moves:  2062 
niter:     2  count:    0  breaks:  0  min_S: 6857.0315  max_S: 6891.3795  S: 6857.0315  ΔS:     -32.3158  moves:   944 
niter:     3  count:    0  breaks:  0  min_S: 6796.0629  max_S: 6891.3795  S: 6796.0629  ΔS:     -60.9685  moves:   740 
niter:     4  count:    0  breaks:  0  min_S: 6796.0629  max_S: 6891.3795  S: 6843.8898  ΔS:      47.8268  moves:   631 
niter:     5  count:    0  breaks:  0  min_S: 6796.0629  max_S: 6891.3795  S: 6854.8522  ΔS:      10.9624  moves:   784 
...
niter:   996  count:    0  breaks:  0  min_S: 6602.7317  max_S: 6912.9081  S: 6798.7139  ΔS:      6.91162  moves:   569 
niter:   997  count:    0  breaks:  0  min_S: 6602.7317  max_S: 6912.9081  S: 6801.4644  ΔS:      2.75049  moves:   929 
niter:   998  count:    0  breaks:  0  min_S: 6602.7317  max_S: 6912.9081  S: 6805.2167  ΔS:      3.75223  moves:  1189 
niter:   999  count:    0  breaks:  0  min_S: 6602.7317  max_S: 6912.9081  S: 6810.0207  ΔS:      4.80406  moves:  1329 
[Mon Jun 16 17:28:57 2025]
Finished job 5.
6 of 9 steps (67%) done
Select jobs to execute...

[Mon Jun 16 17:28:57 2025]
rule exportBlocks:
    input: cache/trimmed_graph/fdr-1e-2/normalized_counts_1.xml.gz, cache/MCMC/blocks/fdr-1e-2/normalized_counts_1.dill
    output: cache/blockSummary/fdr-1e-2/normalized_counts_1/block_summary.csv, cache/blockSummary/fdr-1e-2/normalized_counts_1/gene_
block.csv
    log: logs/exportBlocks/fdr-1e-2/normalized_counts_1.log
    jobid: 2
    reason: Missing output files: cache/blockSummary/fdr-1e-2/normalized_counts_1/block_summary.csv, cache/blockSummary/fdr-1e-2/nor
malized_counts_1/gene_block.csv; Input files updated by another job: cache/trimmed_graph/fdr-1e-2/normalized_counts_1.xml.gz, cache/
MCMC/blocks/fdr-1e-2/normalized_counts_1.dill
    wildcards: fdr=1e-2, label=normalized_counts_1
    threads: 16
    resources: tmpdir=/tmp

/home/chrwhe/micromamba/envs/snakemake_env6/lib/python3.9/site-packages/numpy/_core/fromnumeric.py:3596: RuntimeWarning: Mean of emp
ty slice.
  return _methods._mean(a, axis=axis, dtype=dtype,
/home/chrwhe/micromamba/envs/snakemake_env6/lib/python3.9/site-packages/numpy/_core/_methods.py:138: RuntimeWarning: invalid value e
ncountered in scalar divide
  ret = ret.dtype.type(ret / rcount)
/home/chrwhe/micromamba/envs/snakemake_env6/lib/python3.9/site-packages/numpy/_core/fromnumeric.py:3596: RuntimeWarning: Mean of emp
ty slice.
  return _methods._mean(a, axis=axis, dtype=dtype,
/home/chrwhe/micromamba/envs/snakemake_env6/lib/python3.9/site-packages/numpy/_core/_methods.py:138: RuntimeWarning: invalid value e
ncountered in scalar divide
  ret = ret.dtype.type(ret / rcount)

#####
tree -L 4 cache/                                              
cache/
├── annealedBlock
│   └── fdr-1e-2
│       └── normalized_counts_1.dill
├── blockSummary
│   └── fdr-1e-2
│       └── normalized_counts_1
│           ├── background.csv
│           ├── block_summary.csv
│           ├── gene_block.csv
│           ├── Level_1
│           └── Level_2
├── equilibrate
│   └── fdr-1e-2
│       └── normalized_counts_1.dill
├── graph
│   └── normalized_counts_1.xml.gz
├── initialBlock
│   └── fdr-1e-2
│       └── normalized_counts_1.dill
├── MCMC
│   ├── blocks
│   │   └── fdr-1e-2
│   │       └── normalized_counts_1.dill
│   └── hist
│       └── fdr-1e-2
│           └── normalized_counts_1.dill
└── trimmed_graph
    └── fdr-1e-2
        └── normalized_counts_1.xml.gz
