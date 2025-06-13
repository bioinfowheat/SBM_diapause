

testing_SBM_tools/SBM-tools/snakemake$ tree -L 4 cache/
cache/
├── graph
│   ├── normalized_counts_1_voom.xml.gz
│   ├── normalized_counts_1.xml.gz
│   ├── normalized_counts_2_voom.xml.gz
│   ├── normalized_counts_2.xml.gz
│   ├── normalized_counts_3_voom.xml.gz
│   ├── normalized_counts_3.xml.gz
│   ├── normalized_counts_voom.xml.gz
│   └── normalized_counts.xml.gz
└── trimmed_graph
    └── fdr-1e-4
        ├── normalized_counts_1_voom.xml.gz
        ├── normalized_counts_1.xml.gz
        └── normalized_counts_3.xml.gz




test_run2/SBM-tools/snakemake$ tree -L 4 cache/
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


