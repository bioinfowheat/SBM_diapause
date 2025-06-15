# trying to start a new run with some of the generated files from a previous 
# run. copying those into folders to see if it detects them and starts where they left off

cd /mnt/griffin/chrwhe/SBM_2

micromamba deactivate

# setting environment
mamba create -n snakemake_env6 python=3.9
micromamba activate snakemake_env6
micromamba install -c conda-forge snakemake statsmodels graph-tool scikit-learn multipy dill

critical libmamba filesystem error: cannot create directories: No space left on device 
[/home/chrwhe/micromamba/envs/snakemake_env6/lib/python3.9/site-packages/pysftp-0.2.9-py2.7.egg-info]



mkdir test_run6
cd test_run6

git clone https://github.com/ayroles-lab/SBM-tools.git


cd SBM-tools/rawData/batch
cp /mnt/griffin/chrwhe/SBM_2/Input_SBM_ab/normalized_counts.tsv .
rm example*  
cd ../../snakemake

# existing folder path
mkdir cache
mkdir cache/trimmed_graph/
mkdir cache/trimmed_graph/fdr-1e-2/
mkdir cache/initialBlock
mkdir cache/initialBlock/fdr-1e-2/
mkdir cache/graph/
mkdir cache/annealedBlock
mkdir cache/annealedBlock/fdr-1e-2/
# fill
file_folder=/mnt/griffin/chrwhe/SBM_2/test_run3/SBM-tools/snakemake/cache
cp $file_folder/trimmed_graph/fdr-1e-2/normalized_counts.xml.gz cache/trimmed_graph/fdr-1e-2/
cp $file_folder/initialBlock/fdr-1e-2/normalized_counts.dill cache/initialBlock/fdr-1e-2/
cp $file_folder/graph/normalized_counts.xml.gz cache/graph/
cp $file_folder/annealedBlock/fdr-1e-2/normalized_counts.dill cache/annealedBlock/fdr-1e-2/

tree -L 4 cache
cache
├── annealedBlock
│   └── fdr-1e-2
│       └── normalized_counts.dill
├── graph
│   └── normalized_counts.xml.gz
├── initialBlock
│   └── fdr-1e-2
│       └── normalized_counts.dill
└── trimmed_graph
    └── fdr-1e-2
        └── normalized_counts.xml.gz


nano config.yaml
fdr: [1e-2]
type: "batch" # "layer" or "batch"


snakemake -s Snakefile_batchGraph --cores all &> screen.log

