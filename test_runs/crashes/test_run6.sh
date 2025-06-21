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

# froze at 
tail *log                                                     
niter:   929  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1382006.3  ΔS:      124.698  moves:  7203            
niter:   930  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1381981.5  ΔS:     -24.8116  moves:  2550            
niter:   931  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1382093.2  ΔS:      111.660  moves:  4310            
niter:   932  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1382089.7  ΔS:     -3.54012  moves:  3108            
niter:   933  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1382054.7  ΔS:     -34.9079  moves:  4150            
niter:   934  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1381940.1  ΔS:     -114.692  moves:  5117            
niter:   935  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1381973.6  ΔS:      33.5061  moves:  5134            
niter:   936  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1381872.6  ΔS:     -100.953  moves:  4247            
niter:   937  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1382043.9  ΔS:      171.329  moves:  4799            
niter:   938  count:    0  breaks:  0  min_S: 1381595.0  max_S: 1382960.1  S: 1381941.2  ΔS:     -102.741  moves:  5056 

in the MCMC step.

# trying to get parallel trim to work in this test set.
