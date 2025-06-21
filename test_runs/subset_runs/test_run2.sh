#########
# not sure if having the large files in amongst these smaller test sets will slow things down ... is each step really separate?
/mnt/griffin/chrwhe/SBM_2/test_run2

git clone https://github.com/ayroles-lab/SBM-tools.git

cd SBM-tools/rawData/batch
cp /mnt/griffin/chrwhe/SBM_2/Input_SBM_ab/normalized_*tsv .
rm example* normalized_counts.tsv normalized_counts_voom.tsv
# 6 files left
ls
# normalized_counts_1.tsv       normalized_counts_2.tsv       normalized_counts_3.tsv
# normalized_counts_1_voom.tsv  normalized_counts_2_voom.tsv  normalized_counts_3_voom.tsv

cd ../../snakemake
micromamba activate snakemake_env
snakemake -s Snakefile_batchGraph --cores all &> screen.log

more logs/trim/fdr-1e-4/normalized_counts_2.log
2025-06-11 19:27:10 FDR level:1e-4
2025-06-11 19:27:10 Reading full graph...
2025-06-11 19:27:11 Trimming...
2025-06-11 19:27:11 200 genes
2025-06-11 19:27:11 Density: 0.0
2025-06-11 19:27:11 Uncaught exception: Traceback (most recent call last):
  File "/mnt/griffin/chrwhe/SBM_2/test_run2/SBM-tools/snakemake/.snakemake/scripts/tmpxq80849_.trim_networks.py", line 99, in <modul
e>
    logging.info("Min correlation: " + str(min(gi.edge_properties["spearman"].a)))
ValueError: min() arg is an empty sequence

######
likely due to empty set due to pvalue threshold.

nano config.yaml
fdr: [1e-4]
#fdr: [1e-2,1e-3,1e-4]
type: "batch" # "layer" or "batch"

# change to 
fdr: [1e-2,1e-3]
#fdr: [1e-2,1e-3,1e-4]
type: "batch" # "layer" or "batch"

snakemake -s Snakefile_batchGraph --cores all &> screen2.log

tree -L 4 cache/

cache/
├── annealedBlock
│   ├── fdr-1e-2
│   │   ├── normalized_counts_1.dill
│   │   ├── normalized_counts_1_voom.dill
│   │   ├── normalized_counts_2.dill
│   │   ├── normalized_counts_2_voom.dill
│   │   ├── normalized_counts_3.dill
│   │   └── normalized_counts_3_voom.dill
│   └── fdr-1e-3
│       ├── normalized_counts_1.dill
│       ├── normalized_counts_1_voom.dill
│       ├── normalized_counts_2.dill
│       ├── normalized_counts_2_voom.dill
│       ├── normalized_counts_3.dill
│       └── normalized_counts_3_voom.dill
├── blockSummary
│   └── fdr-1e-2
│       ├── normalized_counts_2
│       │   ├── background.csv
│       │   ├── block_summary.csv
│       │   ├── gene_block.csv
│       │   └── Level_1
│       └── normalized_counts_2_voom
│           ├── background.csv
│           ├── block_summary.csv
│           ├── gene_block.csv
│           └── Level_1
├── equilibrate
│   ├── fdr-1e-2
│   │   ├── normalized_counts_1.dill
│   │   ├── normalized_counts_2.dill
│   │   ├── normalized_counts_2_voom.dill
│   │   ├── normalized_counts_3.dill
│   │   └── normalized_counts_3_voom.dill
│   └── fdr-1e-3
│       ├── normalized_counts_1.dill
│       ├── normalized_counts_1_voom.dill
│       ├── normalized_counts_2.dill
│       ├── normalized_counts_2_voom.dill
│       ├── normalized_counts_3.dill
│       └── normalized_counts_3_voom.dill
├── graph
│   ├── normalized_counts_1_voom.xml.gz
│   ├── normalized_counts_1.xml.gz
│   ├── normalized_counts_2_voom.xml.gz
│   ├── normalized_counts_2.xml.gz
│   ├── normalized_counts_3_voom.xml.gz
│   └── normalized_counts_3.xml.gz
├── initialBlock
│   ├── fdr-1e-2
│   │   ├── normalized_counts_1.dill
│   │   ├── normalized_counts_1_voom.dill
│   │   ├── normalized_counts_2.dill
│   │   ├── normalized_counts_2_voom.dill
│   │   ├── normalized_counts_3.dill
│   │   └── normalized_counts_3_voom.dill
│   ├── fdr-1e-3
│   │   ├── normalized_counts_1.dill
│   │   ├── normalized_counts_1_voom.dill
│   │   ├── normalized_counts_2.dill
│   │   ├── normalized_counts_2_voom.dill
│   │   ├── normalized_counts_3.dill
│   │   └── normalized_counts_3_voom.dill
│   └── fdr-1e-4
│       └── normalized_counts_1.dill
├── MCMC
│   ├── blocks
│   │   ├── fdr-1e-2
│   │   │   ├── normalized_counts_1.dill
│   │   │   ├── normalized_counts_2.dill
│   │   │   ├── normalized_counts_2_voom.dill
│   │   │   └── normalized_counts_3_voom.dill
│   │   └── fdr-1e-3
│   │       └── normalized_counts_1_voom.dill
│   └── hist
│       ├── fdr-1e-2
│       │   ├── normalized_counts_1.dill
│       │   ├── normalized_counts_2.dill
│       │   ├── normalized_counts_2_voom.dill
│       │   └── normalized_counts_3_voom.dill
│       └── fdr-1e-3
│           └── normalized_counts_1_voom.dill
└── trimmed_graph
    ├── fdr-1e-2
    │   ├── normalized_counts_1_voom.xml.gz
    │   ├── normalized_counts_1.xml.gz
    │   ├── normalized_counts_2_voom.xml.gz
    │   ├── normalized_counts_2.xml.gz
    │   ├── normalized_counts_3_voom.xml.gz
    │   └── normalized_counts_3.xml.gz
    ├── fdr-1e-3
    │   ├── normalized_counts_1_voom.xml.gz
    │   ├── normalized_counts_1.xml.gz
    │   ├── normalized_counts_2_voom.xml.gz
    │   ├── normalized_counts_2.xml.gz
    │   ├── normalized_counts_3_voom.xml.gz
    │   └── normalized_counts_3.xml.gz
    └── fdr-1e-4
        ├── normalized_counts_1_voom.xml.gz
        ├── normalized_counts_1.xml.gz
        └── normalized_counts_3.xml.gz

tree -L 4 logs/
tree -L 4 logs/
logs/
├── annealedBlock
│   ├── fdr-1e-2
│   │   ├── normalized_counts_1.log
│   │   ├── normalized_counts_1_voom.log
│   │   ├── normalized_counts_2.log
│   │   ├── normalized_counts_2_voom.log
│   │   ├── normalized_counts_3.log
│   │   └── normalized_counts_3_voom.log
│   └── fdr-1e-3
│       ├── normalized_counts_1.log
│       ├── normalized_counts_1_voom.log
│       ├── normalized_counts_2.log
│       ├── normalized_counts_2_voom.log
│       ├── normalized_counts_3.log
│       └── normalized_counts_3_voom.log
├── equilibrate
│   ├── fdr-1e-2
│   │   ├── normalized_counts_1.log
│   │   ├── normalized_counts_2.log
│   │   ├── normalized_counts_2_voom.log
│   │   ├── normalized_counts_3.log
│   │   └── normalized_counts_3_voom.log
│   └── fdr-1e-3
│       ├── normalized_counts_1.log
│       ├── normalized_counts_1_voom.log
│       ├── normalized_counts_2.log
│       ├── normalized_counts_2_voom.log
│       ├── normalized_counts_3.log
│       └── normalized_counts_3_voom.log
├── exportBlocks
│   └── fdr-1e-2
│       ├── normalized_counts_2.log
│       └── normalized_counts_2_voom.log
├── GO
│   └── fdr-1e-2
│       ├── normalized_counts_2.log
│       └── normalized_counts_2_voom.log
├── intialBlock
│   ├── fdr-1e-2
│   │   ├── normalized_counts_1.log
│   │   ├── normalized_counts_1_voom.log
│   │   ├── normalized_counts_2.log
│   │   ├── normalized_counts_2_voom.log
│   │   ├── normalized_counts_3.log
│   │   └── normalized_counts_3_voom.log
│   ├── fdr-1e-3
│   │   ├── normalized_counts_1.log
│   │   ├── normalized_counts_1_voom.log
│   │   ├── normalized_counts_2.log
│   │   ├── normalized_counts_2_voom.log
│   │   ├── normalized_counts_3.log
│   │   └── normalized_counts_3_voom.log
│   └── fdr-1e-4
│       └── normalized_counts_1.log
├── MCMC
│   ├── fdr-1e-2
│   │   ├── normalized_counts_1.log
│   │   ├── normalized_counts_2.log
│   │   ├── normalized_counts_2_voom.log
│   │   └── normalized_counts_3_voom.log
│   └── fdr-1e-3
│       └── normalized_counts_1_voom.log
├── slurm
│   └── README.md
└── trim
    ├── fdr-1e-2
    │   ├── normalized_counts_1.log
    │   ├── normalized_counts_1_voom.log
    │   ├── normalized_counts_2.log
    │   ├── normalized_counts_2_voom.log
    │   ├── normalized_counts_3.log
    │   └── normalized_counts_3_voom.log
    ├── fdr-1e-3
    │   ├── normalized_counts_1.log
    │   ├── normalized_counts_1_voom.log
    │   ├── normalized_counts_2.log
    │   ├── normalized_counts_2_voom.log
    │   ├── normalized_counts_3.log
    │   └── normalized_counts_3_voom.log
    └── fdr-1e-4
        ├── normalized_counts_1.log
        ├── normalized_counts_1_voom.log
        ├── normalized_counts_2.log
        └── normalized_counts_3.log


# of the 6 dataset, and three fdr levels: [1e-2,1e-3,1e-4]
# looks like only 1 dataset at the least restrictive pvalue
# came through normalized_counts_2

zip  

######
# collect and zip necessary files

rawdatafile=normalized_counts_2
zip run2.$rawdatafile.zip \
cache/blockSummary/fdr-1e-2/$rawdatafile/gene_block.csv \
cache/blockSummary/fdr-1e-2/$rawdatafile/block_summary.csv \
cache/trimmed_graph/fdr-1e-2/$rawdatafile.xml.gz \
cache/MCMC/blocks/fdr-1e-2/$rawdatafile.dill

rawdatafile=normalized_counts_2_voom
zip run2.$rawdatafile.zip \
cache/blockSummary/fdr-1e-2/$rawdatafile/gene_block.csv \
cache/blockSummary/fdr-1e-2/$rawdatafile/block_summary.csv \
cache/trimmed_graph/fdr-1e-2/$rawdatafile.xml.gz \
cache/MCMC/blocks/fdr-1e-2/$rawdatafile.dill


scp chrwhe@duke.zoologi.su.se:/mnt/griffin/chrwhe/SBM_2/test_run2/SBM-tools/snakemake/run2.\*.zip .

