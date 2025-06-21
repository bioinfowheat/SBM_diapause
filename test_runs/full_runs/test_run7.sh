# trying to start a new run with some of the generated files from a previous 
# run. copying those into folders to see if it detects them and starts where they left off

cd /mnt/griffin/chrwhe/SBM_2

micromamba deactivate

# setting environment
mamba create -n snakemake_env7 python=3.9
micromamba activate snakemake_env7
micromamba install -c conda-forge snakemake statsmodels graph-tool scikit-learn multipy dill



mkdir test_run7
cd test_run7

git clone https://github.com/ayroles-lab/SBM-tools.git


cd SBM-tools/rawData/batch
cp /mnt/griffin/chrwhe/SBM_2/Input_SBM_ab/normalized_counts.tsv .
rm example*  
cd ../../snakemake

# existing folder path
mkdir cache
mkdir cache/trimmed_graph/
mkdir cache/trimmed_graph/fdr-1e-3/
mkdir cache/initialBlock
mkdir cache/initialBlock/fdr-1e-3/
mkdir cache/graph/
mkdir cache/annealedBlock
mkdir cache/annealedBlock/fdr-1e-3/
# fill
file_folder=/mnt/griffin/chrwhe/SBM_2/test_run3/SBM-tools/snakemake/cache
cp $file_folder/trimmed_graph/fdr-1e-3/normalized_counts.xml.gz cache/trimmed_graph/fdr-1e-3/
cp $file_folder/initialBlock/fdr-1e-3/normalized_counts.dill cache/initialBlock/fdr-1e-3/
cp $file_folder/graph/normalized_counts.xml.gz cache/graph/
cp $file_folder/annealedBlock/fdr-1e-3/normalized_counts.dill cache/annealedBlock/fdr-1e-3/


nano config.yaml
fdr: [1e-3]
type: "batch" # "layer" or "batch"


snakemake -s Snakefile_batchGraph --cores all &> screen.log

has been stuck at 
test_run7/SBM-tools/snakemake$ tail *log
niter:   929  count:    0  breaks:  0  min_S: 150570.67  max_S: 155029.58  S: 152947.58  ΔS:      29.5657  moves: 22919 
niter:   930  count:    0  breaks:  0  min_S: 150570.67  max_S: 155029.58  S: 152820.62  ΔS:     -126.956  moves: 25587 
niter:   931  count:    0  breaks:  0  min_S: 150570.67  max_S: 155029.58  S: 152756.61  ΔS:     -64.0128  moves: 21215 
niter:   932  count:    0  breaks:  0  min_S: 150570.67  max_S: 155029.58  S: 152850.88  ΔS:      94.2695  moves: 21444 
niter:   933  count:    0  breaks:  0  min_S: 150570.67  max_S: 155029.58  S: 152681.26  ΔS:     -169.614  moves: 23294 
niter:   934  count:    0  breaks:  0  min_S: 150570.67  max_S: 155029.58  S: 152602.71  ΔS:     -78.5530  moves: 21290 
niter:   935  count:    0  breaks:  0  min_S: 150570.67  max_S: 155029.58  S: 152561.66  ΔS:     -41.0500  moves: 21205 
niter:   936  count:    0  breaks:  0  min_S: 150570.67  max_S: 155029.58  S: 152487.73  ΔS:     -73.9280  moves: 24059 
niter:   937  count:    0  breaks:  0  min_S: 150570.67  max_S: 155029.58  S: 152214.75  ΔS:     -272.986  moves: 20194 
niter:   938  count:    0  breaks:  0  min_S: 150570.67  max_S: 155029.58  S: 152400.62  ΔS:      185.870  moves: 24029 

for 3 days. Killing.


######
# restarting

run_folder=test_run7
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

cp /mnt/griffin/chrwhe/SBM_2/SBM-tools/snakemake/scripts/exportBlocks.py /mnt/griffin/chrwhe/SBM_2/test_run7/SBM-tools/snakemake/scripts


# sanity check 
more /mnt/griffin/chrwhe/SBM_2/$run_folder/SBM-tools/snakemake/scripts/make_stacked_networks.py

# setting environment
micromamba activate snakemake_env7
# nano config to 
nano config.yaml
fdr: [1e-3]

# try runnning.
snakemake -s Snakefile_batchGraph --cores all &> screen.log
