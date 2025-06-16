
# go to the folder you want to set up a test run
cd /mnt/griffin/XXX/testruns

# setting environment
mamba create -n snakemake_env1 python=3.9
micromamba activate snakemake_env1
micromamba install -c conda-forge snakemake statsmodels graph-tool scikit-learn multipy dill

# if you see the error below we have run out of space on the home folder and need to clean things
# email me ... 
# critical libmamba filesystem error: cannot create directories: No space left on device 
# [/home/chrwhe/micromamba/envs/snakemake_env6/lib/python3.9/site-packages/pysftp-0.2.9-py2.7.egg-info]



mkdir test_run1
cd test_run1

git clone https://github.com/ayroles-lab/SBM-tools.git


cd SBM-tools/rawData/batch
# put your dataset into this folder
cp /mnt/griffin/chrwhe/SBM_2/Input_SBM500/Input_SBM_ab500/normalized_counts_1.tsv .
# remove the example file
rm example*  
# go back to run folder
cd ../../snakemake

# nano config file to run at the FDR settting you want
# default is 1e-4
nano config.yaml
fdr: [1e-2]

# try runnning. This will dump everything that would normally go to the screen into this "screen.log" file (could be any name you want)
snakemake -s Snakefile_batchGraph --cores all &> screen.log

