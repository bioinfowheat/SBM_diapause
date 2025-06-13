######################
# test_run3
######################
rawData/batch$ ls
normalized_counts.tsv  normalized_counts_voom.tsv


head config.yaml 
fdr: [1e-2,1e-3,1e-4]
type: "batch" # "layer" or "batch"


test_run3/SBM-tools/snakemake$ tree -L 4 cache/
cache/
├── annealedBlock
│   ├── fdr-1e-2
│   │   └── normalized_counts.dill
│   ├── fdr-1e-3
│   │   └── normalized_counts.dill
│   └── fdr-1e-4
│       └── normalized_counts.dill
├── blockSummary
│   └── fdr-1e-4
│       └── normalized_counts
│           ├── background.csv
│           └── Level_1
├── equilibrate
│   ├── fdr-1e-3
│   │   └── normalized_counts.dill
│   └── fdr-1e-4
│       └── normalized_counts.dill
├── graph
│   └── normalized_counts.xml.gz
├── initialBlock
│   ├── fdr-1e-2
│   │   └── normalized_counts.dill
│   ├── fdr-1e-3
│   │   └── normalized_counts.dill
│   └── fdr-1e-4
│       └── normalized_counts.dill
├── MCMC
│   ├── blocks
│   │   └── fdr-1e-4
│   │       └── normalized_counts.dill
│   └── hist
│       └── fdr-1e-4
│           └── normalized_counts.dill
└── trimmed_graph
    ├── fdr-1e-2
    │   └── normalized_counts.xml.gz
    ├── fdr-1e-3
    │   └── normalized_counts.xml.gz
    └── fdr-1e-4
        └── normalized_counts.xml.gz


Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 64
Rules claiming more threads will be scaled down.
Job stats:
job             count
------------  -------
GO                  6
MCMC                6
all                 1
annealing           6
equilibrate         6
exportBlocks        6
graph               2
minSBM              6
trim                6
total              45

Select jobs to execute...
Failed to solve scheduling problem with ILP solver. Falling back to greedy solver. Run Snakemake with --verbose to see the full solv
er output for debugging the problem.

[Wed Jun 11 23:43:32 2025]
rule graph:
    input: ../rawData/batch/normalized_counts_voom.tsv
    output: cache/graph/normalized_counts_voom.xml.gz
    jobid: 12
    reason: Missing output files: cache/graph/normalized_counts_voom.xml.gz
    wildcards: label=normalized_counts_voom
    resources: tmpdir=/tmp

[Wed Jun 11 23:43:32 2025]
rule graph:
    input: ../rawData/batch/normalized_counts.tsv
    output: cache/graph/normalized_counts.xml.gz
    jobid: 4
    reason: Missing output files: cache/graph/normalized_counts.xml.gz
    wildcards: label=normalized_counts
    resources: tmpdir=/tmp


/home/chrwhe/micromamba/envs/snakemake_env/lib/python3.9/site-packages/numpy/core/_methods.py:189: RuntimeWarning: invalid value enc
ountered in double_scalars
  ret = ret.dtype.type(ret / rcount)
/home/chrwhe/micromamba/envs/snakemake_env/lib/python3.9/site-packages/numpy/core/fromnumeric.py:3474: RuntimeWarning: Mean of empty
 slice.
  return _methods._mean(a, axis=axis, dtype=dtype,

/home/chrwhe/micromamba/envs/snakemake_env/lib/python3.9/site-packages/numpy/core/fromnumeric.py:3474: RuntimeWarning: Mean of empty
 slice.
  return _methods._mean(a, axis=axis, dtype=dtype,
/home/chrwhe/micromamba/envs/snakemake_env/lib/python3.9/site-packages/numpy/core/_methods.py:189: RuntimeWarning: invalid value enc
ountered in double_scalars
  ret = ret.dtype.type(ret / rcount)

[Thu Jun 12 22:06:19 2025]
Error in rule exportBlocks:
    jobid: 32
    input: cache/trimmed_graph/fdr-1e-4/normalized_counts.xml.gz, cache/MCMC/blocks/fdr-1e-4/normalized_counts.dill
    output: cache/blockSummary/fdr-1e-4/normalized_counts/block_summary.csv, cache/blockSummary/fdr-1e-4/normalized_counts/gene_bloc
k.csv
    log: logs/exportBlocks/fdr-1e-4/normalized_counts.log (check log file(s) for error details)

RuleException:
CalledProcessError in file /mnt/griffin/chrwhe/SBM_2/test_run3/SBM-tools/snakemake/Snakefile_batchGraph, line 131:
Command 'set -euo pipefail;  /home/chrwhe/micromamba/envs/snakemake_env/bin/python3.9 /mnt/griffin/chrwhe/SBM_2/test_run3/SBM-tools/
snakemake/.snakemake/scripts/tmpuuzrg7et.exportBlocks.py' returned non-zero exit status 1.
  File "/mnt/griffin/chrwhe/SBM_2/test_run3/SBM-tools/snakemake/Snakefile_batchGraph", line 131, in __rule_exportBlocks
  File "/home/chrwhe/micromamba/envs/snakemake_env/lib/python3.9/concurrent/futures/thread.py", line 58, in run
[Thu Jun 12 22:49:50 2025]

more logs/exportBlocks/fdr-1e-4/normalized_counts.log
2025-06-12 22:05:35 Reading trimmed graph...
2025-06-12 22:05:36 Arctan transform on correlations...
2025-06-12 22:05:36 Loading blocks...
2025-06-12 22:05:36 Creating nested block model...
2025-06-12 22:05:40 Creating block DataFrame...
2025-06-12 22:06:18 Calculating block sizes...
2025-06-12 22:06:19 Creating gene lists...
2025-06-12 22:06:19 At level: 1
2025-06-12 22:06:19 Uncaught exception: Traceback (most recent call last):
  File "/mnt/griffin/chrwhe/SBM_2/test_run3/SBM-tools/snakemake/.snakemake/scripts/tmpuuzrg7et.exportBlocks.py", line 131, in <modul
e>
    sys.error("q_r is larger than one.")
AttributeError: module 'sys' has no attribute 'error'

tree -L 4 logs/
logs/
├── annealedBlock
│   ├── fdr-1e-2
│   │   └── normalized_counts.log
│   ├── fdr-1e-3
│   │   └── normalized_counts.log
│   └── fdr-1e-4
│       └── normalized_counts.log
├── equilibrate
│   ├── fdr-1e-3
│   │   └── normalized_counts.log
│   └── fdr-1e-4
│       └── normalized_counts.log
├── exportBlocks
│   └── fdr-1e-4
│       └── normalized_counts.log
├── intialBlock
│   ├── fdr-1e-2
│   │   └── normalized_counts.log
│   ├── fdr-1e-3
│   │   └── normalized_counts.log
│   └── fdr-1e-4
│       └── normalized_counts.log
├── MCMC
│   └── fdr-1e-4
│       └── normalized_counts.log
├── slurm
│   └── README.md
└── trim
    ├── fdr-1e-2
    │   └── normalized_counts.log
    ├── fdr-1e-3
    │   └── normalized_counts.log
    └── fdr-1e-4
        └── normalized_counts.log


# this is still running ... 

# going to start another instance.

cp -r SBM-tools SBM-tools_rerun
cd /mnt/griffin/chrwhe/SBM_2/test_run3/SBM-tools_rerun
# removed voom data
snakemake -s Snakefile_batchGraph --cores all &> screen_rerun.log