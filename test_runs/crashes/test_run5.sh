# worried about invoking a single environement twice, so I'm making a new one here

cd /mnt/griffin/chrwhe/SBM_2

micromamba deactivate

# setting environment
mamba create -n snakemake_env3 python=3.9
micromamba activate snakemake_env3
micromamba install -c conda-forge snakemake statsmodels graph-tool scikit-learn multipy dill

mkdir test_run5
cd test_run5

git clone https://github.com/ayroles-lab/SBM-tools.git


cd SBM-tools/rawData/batch
cp /mnt/griffin/chrwhe/SBM_2/Input_SBM_ab/normalized_counts_voom.tsv .
rm example*  
cd ../../snakemake

nano config.yaml
fdr: [1e-2,1e-3,1e-4]
type: "batch" # "layer" or "batch"


snakemake -s Snakefile_batchGraph --cores all &> screen.log

more /mnt/griffin/chrwhe/SBM_2/test_run5/SBM-tools/snakemake/screen.log
Job stats:
job             count
------------  -------
GO                  3
MCMC                3
all                 1
annealing           3
equilibrate         3
exportBlocks        3
graph               1
minSBM              3
trim                3
total              23

Select jobs to execute...

[Fri Jun 13 14:41:01 2025]
rule graph:
    input: ../rawData/batch/normalized_counts.tsv
    output: cache/graph/normalized_counts.xml.gz
    jobid: 4
    reason: Missing output files: cache/graph/normalized_counts.xml.gz
    wildcards: label=normalized_counts
    resources: tmpdir=/tmp
    
tree -L 4 cache/
cache/
├── annealedBlock
│   ├── fdr-1e-2
│   │   └── normalized_counts_voom.dill
│   ├── fdr-1e-3
│   │   └── normalized_counts_voom.dill
│   └── fdr-1e-4
│       └── normalized_counts_voom.dill
├── blockSummary
│   └── fdr-1e-4
│       └── normalized_counts_voom
│           ├── background.csv
│           └── Level_1
├── equilibrate
│   ├── fdr-1e-3
│   │   └── normalized_counts_voom.dill
│   └── fdr-1e-4
│       └── normalized_counts_voom.dill
├── graph
│   └── normalized_counts_voom.xml.gz
├── initialBlock
│   ├── fdr-1e-2
│   │   └── normalized_counts_voom.dill
│   ├── fdr-1e-3
│   │   └── normalized_counts_voom.dill
│   └── fdr-1e-4
│       └── normalized_counts_voom.dill
├── MCMC
│   ├── blocks
│   │   └── fdr-1e-4
│   │       └── normalized_counts_voom.dill
│   └── hist
│       └── fdr-1e-4
│           └── normalized_counts_voom.dill
└── trimmed_graph
    ├── fdr-1e-2
    │   └── normalized_counts_voom.xml.gz
    ├── fdr-1e-3
    │   └── normalized_counts_voom.xml.gz
    └── fdr-1e-4
        └── normalized_counts_voom.xml.gz