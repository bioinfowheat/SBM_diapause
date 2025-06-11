
# setting environment
mamba create -n snakemake_env python=3.9
micromamba activate snakemake_env
micromamba install -c conda-forge snakemake statsmodels graph-tool scikit-learn multipy dill


git clone https://github.com/ayroles-lab/SBM-tools.git

cd rawData/batch
cp ../../../../Input_SBM_ab/normalized_counts_1.tsv ../../../../Input_SBM_ab/normalized_counts_2.tsv .
#
ls
# example.tsv  normalized_counts_1.tsv  normalized_counts_2.tsv

cd /mnt/griffin/chrwhe/SBM_2/testing_SBM_tools/SBM-tools/snakemake
micromamba activate snakemake_env
snakemake -s Snakefile_batchGraph --cores all


##### 
# wow, it runs all datasets in the rawData/batch folder in parallel .... VERY NICE

cd /mnt/griffin/chrwhe/SBM_2/testing_SBM_tools/SBM-tools/rawData/batch
cp ../../../../Input_SBM_ab/n*tsv .
rm example*

cd /mnt/griffin/chrwhe/SBM_2/testing_SBM_tools/SBM-tools/snakemake
snakemake -s Snakefile_batchGraph --cores all

Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 64
Rules claiming more threads will be scaled down.
Job stats:
job             count
------------  -------
GO                  8
MCMC                8
all                 1
annealing           8
equilibrate         8
exportBlocks        8
graph               8
minSBM              8
trim                8
total              65

Select jobs to execute...
Failed to solve scheduling problem with ILP solver. Falling back to greedy solver. Run Snakemake with --verbose to see the full solver output for debugging the problem.

[Wed Jun 11 19:15:58 2025]
rule graph:
    input: ../rawData/batch/normalized_counts_2_voom.tsv
    output: cache/graph/normalized_counts_2_voom.xml.gz
    jobid: 28
    reason: Missing output files: cache/graph/normalized_counts_2_voom.xml.gz
    wildcards: label=normalized_counts_2_voom
    resources: tmpdir=/tmp

[Wed Jun 11 19:15:58 2025]
rule graph:
    input: ../rawData/batch/normalized_counts_1.tsv
    output: cache/graph/normalized_counts_1.xml.gz
    jobid: 4
    reason: Missing output files: cache/graph/normalized_counts_1.xml.gz
    wildcards: label=normalized_counts_1
    resources: tmpdir=/tmp

[Wed Jun 11 19:15:58 2025]
rule graph:
    input: ../rawData/batch/normalized_counts_3_voom.tsv
    output: cache/graph/normalized_counts_3_voom.xml.gz
    jobid: 44
    reason: Missing output files: cache/graph/normalized_counts_3_voom.xml.gz
    wildcards: label=normalized_counts_3_voom
    resources: tmpdir=/tmp

[Wed Jun 11 19:15:58 2025]
rule graph:
    input: ../rawData/batch/normalized_counts_voom.tsv
    output: cache/graph/normalized_counts_voom.xml.gz
    jobid: 60
    reason: Missing output files: cache/graph/normalized_counts_voom.xml.gz
    wildcards: label=normalized_counts_voom
    resources: tmpdir=/tmp

[Wed Jun 11 19:15:58 2025]
rule graph:
    input: ../rawData/batch/normalized_counts_1_voom.tsv
    output: cache/graph/normalized_counts_1_voom.xml.gz
    jobid: 20
    reason: Missing output files: cache/graph/normalized_counts_1_voom.xml.gz
    wildcards: label=normalized_counts_1_voom
    resources: tmpdir=/tmp

[Wed Jun 11 19:15:58 2025]
rule graph:
    input: ../rawData/batch/normalized_counts_3.tsv
    output: cache/graph/normalized_counts_3.xml.gz
    jobid: 36
    reason: Missing output files: cache/graph/normalized_counts_3.xml.gz
    wildcards: label=normalized_counts_3
    resources: tmpdir=/tmp

[Wed Jun 11 19:15:58 2025]
rule graph:
    input: ../rawData/batch/normalized_counts.tsv
    output: cache/graph/normalized_counts.xml.gz
    jobid: 52
    reason: Missing output files: cache/graph/normalized_counts.xml.gz
    wildcards: label=normalized_counts
    resources: tmpdir=/tmp

[Wed Jun 11 19:15:59 2025]
rule graph:
    input: ../rawData/batch/normalized_counts_2.tsv
    output: cache/graph/normalized_counts_2.xml.gz
    jobid: 12
    reason: Missing output files: cache/graph/normalized_counts_2.xml.gz
    wildcards: label=normalized_counts_2
    resources: tmpdir=/tmp

[Wed Jun 11 19:16:22 2025]
Finished job 36.
1 of 65 steps (2%) done
Select jobs to execute...

[Wed Jun 11 19:16:22 2025]
rule trim:
    input: cache/graph/normalized_counts_3.xml.gz
    output: cache/trimmed_graph/fdr-1e-4/normalized_counts_3.xml.gz
    log: logs/trim/fdr-1e-4/normalized_counts_3.log
    jobid: 35
    reason: Missing output files: cache/trimmed_graph/fdr-1e-4/normalized_counts_3.xml.gz; Input files updated by another job: cache/graph/normalized_counts_3.xml.gz
    wildcards: fdr=1e-4, label=normalized_counts_3
    resources: tmpdir=/tmp

[Wed Jun 11 19:16:22 2025]
Finished job 12.
2 of 65 steps (3%) done
Select jobs to execute...

[Wed Jun 11 19:16:22 2025]
rule trim:
    input: cache/graph/normalized_counts_2.xml.gz
    output: cache/trimmed_graph/fdr-1e-4/normalized_counts_2.xml.gz
    log: logs/trim/fdr-1e-4/normalized_counts_2.log
    jobid: 11
    reason: Missing output files: cache/trimmed_graph/fdr-1e-4/normalized_counts_2.xml.gz; Input files updated by another job: cache/graph/normalized_counts_2.xml.gz
    wildcards: fdr=1e-4, label=normalized_counts_2
    resources: tmpdir=/tmp

[Wed Jun 11 19:16:24 2025]
Finished job 20.
3 of 65 steps (5%) done
Select jobs to execute...

[Wed Jun 11 19:16:24 2025]
rule trim:
    input: cache/graph/normalized_counts_1_voom.xml.gz
    output: cache/trimmed_graph/fdr-1e-4/normalized_counts_1_voom.xml.gz
    log: logs/trim/fdr-1e-4/normalized_counts_1_voom.log
    jobid: 19
    reason: Missing output files: cache/trimmed_graph/fdr-1e-4/normalized_counts_1_voom.xml.gz; Input files updated by another job: cache/graph/normalized_counts_1_voom.xml.gz
    wildcards: fdr=1e-4, label=normalized_counts_1_voom
    resources: tmpdir=/tmp

[Wed Jun 11 19:16:25 2025]
Finished job 4.
4 of 65 steps (6%) done
Select jobs to execute...

[Wed Jun 11 19:16:25 2025]
rule trim:
    input: cache/graph/normalized_counts_1.xml.gz
    output: cache/trimmed_graph/fdr-1e-4/normalized_counts_1.xml.gz
    log: logs/trim/fdr-1e-4/normalized_counts_1.log
    jobid: 3
    reason: Missing output files: cache/trimmed_graph/fdr-1e-4/normalized_counts_1.xml.gz; Input files updated by another job: cache/graph/normalized_counts_1.xml.gz
    wildcards: fdr=1e-4, label=normalized_counts_1
    resources: tmpdir=/tmp

    Error in rule trim:
    jobid: 11
    input: cache/graph/normalized_counts_2.xml.gz
    output: cache/trimmed_graph/fdr-1e-4/normalized_counts_2.xml.gz
    log: logs/trim/fdr-1e-4/normalized_counts_2.log (check log file(s) for error details)

RuleException:
CalledProcessError in file /mnt/griffin/chrwhe/SBM_2/testing_SBM_tools/SBM-tools/snakemake/Snakefile_batchGraph, line 46:
Command 'set -euo pipefail;  /home/chrwhe/micromamba/envs/snakemake_env/bin/python3.9 /mnt/griffin/chrwhe/SBM_2/testing_SBM_tools/SBM-tools/snakemake/.snakemake/scripts/tmp2r17mklv.trim_networks.py' returned non-zero exit status 1.
  File "/mnt/griffin/chrwhe/SBM_2/testing_SBM_tools/SBM-tools/snakemake/Snakefile_batchGraph", line 46, in __rule_trim
  File "/home/chrwhe/micromamba/envs/snakemake_env/lib/python3.9/concurrent/futures/thread.py", line 58, in run


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

