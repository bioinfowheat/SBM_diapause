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

