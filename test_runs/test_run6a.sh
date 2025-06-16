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


