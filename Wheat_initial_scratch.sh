# assessing SBM github

# https://github.com/ayroles-lab/SBM-tools/tree/main

# made a local clone of this
git clone https://github.com/ayroles-lab/SBM-tools.git



# following
# https://github.com/deto/Snakemake_Tutorial

# pip install snakemake
# WARNING: The scripts snakemake and snakemake-bash-completion are installed in
# '/Users/chriswheat/Library/Python/3.9/bin' which is not on PATH.

WARNING: The script tabulate is installed in '/Users/chriswheat/Library/Python/3.9/bin' which is not on PATH.
 Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
  ━━━━━━━━━━━━━━━╸━━━━━━━━━━━━━━━━━━━━━━━━ 14/36 [pulp]  WARNING: The script pulptest is installed in '/Users/chriswheat/Library/Python/3.9/bin' which is not on PATH.
 Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
  ━━━━━━━━━━━━━━━━╸━━━━━━━━━━━━━━━━━━━━━━━ 15/36 [psutil]  WARNING: The script humanfriendly is installed in '/Users/chriswheat/Library/Python/3.9/bin' which is not on PATH.
 Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
  ━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━ 21/36 [docutils]  WARNING: The scripts docutils, rst2html, rst2html4, rst2html5, rst2latex, rst2man, rst2odt, rst2pseudoxml, rst2s5, rst2xetex and rst2xml are installed in '/Users/chriswheat/Library/Python/3.9/bin' which is not on PATH.
 Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
 WARNING: The script yte is installed in '/Users/chriswheat/Library/Python/3.9/bin' which is not on PATH.
 Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━╸━━━━━━━━━━━ 26/36 [smart_open]  WARNING: The scripts jupyter, jupyter-migrate and jupyter-troubleshoot are installed in '/Users/chriswheat/Library/Python/3.9/bin' which is not on PATH.
 Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╸━━━━ 32/36 [gitpython]  WARNING: The script jsonschema is installed in '/Users/chriswheat/Library/Python/3.9/bin' which is not on PATH.
 Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
 WARNING: The script jupyter-trust is installed in '/Users/chriswheat/Library/Python/3.9/bin' which is not on PATH.
 Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╸━ 35/36 [snakemake]  WARNING: The scripts snakemake and snakemake-bash-completion are installed in '/Users/chriswheat/Library/Python/3.9/bin' which is not on PATH.

export PATH=/Users/chriswheat/Library/Python/3.9/bin/:$PATH
export PATH=/Users/chriswheat/Library/Python/3.9/lib/python/site-packages/:$PATH

snakemake Snakefile_batchGraph
###### the above did not work on my local computer

# on server
cd software
# made a local clone of this
git clone https://github.com/ayroles-lab/SBM-tools.git

mamba create -n snakemake_env python=3.9
micromamba activate snakemake_env
micromamba install -c conda-forge snakemake statsmodels graph-tool scikit-learn multipy dill

# note that the graph-tool installed here is MP enabled.

>>> import graph_tool.all as gt
>>> import pandas as pd
>>> import numpy as np
>>> import dill
>>> gt.openmp_enabled()
True



sudo R
	install.packages("log4r", dependencies = TRUE)
	BiocManager::install("org.Hs.eg.db")
	BiocManager::install("clusterProfiler")
# install.packages("clusterProfiler", dependencies = TRUE)
# devtools::install_version("clusterProfiler", version = "4.0.5")
BiocManager::install(version = "3.13")  # Compatible with R 4.1.0
BiocManager::install("clusterProfiler")

#####
> BiocManager::install(version = "3.13")
'getOption("repos")' replaces Bioconductor standard repositories, see
'help("repositories", package = "BiocManager")' for details.
Replacement repositories:
    CRAN: https://cran.mirror.garr.it/CRAN
Warning: unable to access index for repository https://cran.mirror.garr.it/CRAN/src/contrib:
  cannot open URL 'https://cran.mirror.garr.it/CRAN/src/contrib/PACKAGES'
Downgrade 31 packages to Bioconductor version '3.13'? [y/n]: y
Bioconductor version 3.13 (BiocManager 1.30.26), R 4.1.0 (2021-05-18)
Warning messages:
1: package(s) not installed when version(s) same as or greater than current; use
  `force = TRUE` to re-install: 'BiocVersio

# installing dependencies
BiocManager::install(c("DOSE", "enrichplot", "GOSemSim", "GO.db", "AnnotationDbi"))

# Download and install from archived version
url <- "https://bioconductor.org/packages/3.13/bioc/src/contrib/clusterProfiler_4.0.5.tar.gz"
install.packages(url, repos = NULL, type = "source")

install.packages("ggtree", dependencies = TRUE)

#######

snakemake -s Snakefile_batchGraph --cores all
['example.tsv']
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 30
Rules claiming more threads will be scaled down.
Job stats:
job             count
------------  -------
GO                  1
MCMC                1
all                 1
annealing           1
equilibrate         1
exportBlocks        1
graph               1
minSBM              1
trim                1
total               9

Select jobs to execute...

[Sun Jun  8 17:49:39 2025]
rule graph:
    input: ../rawData/batch/example.tsv
    output: cache/graph/example.xml.gz
    jobid: 4
    reason: Missing output files: cache/graph/example.xml.gz
    wildcards: label=example
    resources: tmpdir=/tmp

[Sun Jun  8 17:59:03 2025]
Finished job 4.
1 of 9 steps (11%) done
Select jobs to execute...

[Sun Jun  8 17:59:03 2025]
rule trim:
    input: cache/graph/example.xml.gz
    output: cache/trimmed_graph/fdr-1e-4/example.xml.gz
    log: logs/trim/fdr-1e-4/example.log
    jobid: 3
    reason: Missing output files: cache/trimmed_graph/fdr-1e-4/example.xml.gz; Input files updated by another job: cache/graph/example.xml.gz
    wildcards: fdr=1e-4, label=example
    resources: tmpdir=/tmp

[Sun Jun  8 17:59:05 2025]
Error in rule trim:
    jobid: 3
    input: cache/graph/example.xml.gz
    output: cache/trimmed_graph/fdr-1e-4/example.xml.gz
    log: logs/trim/fdr-1e-4/example.log (check log file(s) for error details)

RuleException:
CalledProcessError in file /mnt/griffin/chrwhe/software/SBM-tools/snakemake/Snakefile_batchGraph, line 46:
Command 'set -euo pipefail;  /home/chrwhe/micromamba/envs/snakemake_env/bin/python3.9 /mnt/griffin/chrwhe/software/SBM-tools/snakemake/.snakemake/scripts/tmpttoy5gf9.trim_networks.py' returned non-zero exit status 1.
  File "/mnt/griffin/chrwhe/software/SBM-tools/snakemake/Snakefile_batchGraph", line 46, in __rule_trim
  File "/home/chrwhe/micromamba/envs/snakemake_env/lib/python3.9/concurrent/futures/thread.py", line 58, in run
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: .snakemake/log/2025-06-08T174938.849667.snakemake.log

# troubleshooting ...
tail logs/trim/fdr-1e-4/example.log
# ModuleNotFoundError: No module named 'sklearn'
# claude says install -c conda-forge scikit-learn

# keep hitting missing components
tail logs/intialBlock/fdr-1e-4/example.log
# ModuleNotFoundError: No module named 'dill'

# restarting and things are runnign well.
snakemake -s Snakefile_batchGraph --cores all

# currently crashes on GO analysis
  there is no package called ‘clusterProfiler’

# once done, which means that all of the required modules have installed, I
# should clean the written files and run clean, then try layered program.

# clean run "
cd /mnt/griffin/chrwhe/SBM_2
# made a local clone of this
git clone https://github.com/ayroles-lab/SBM-tools.git
micromamba activate snakemake_env
cd SBM-tools/snakemake
nano config.yaml
fdr: [1e-5]
snakemake -s Snakefile_batchGraph --cores all

[Sun Jun  8 22:11:20 2025]
Finished job 3.
2 of 9 steps (22%) done
Select jobs to execute...

[Sun Jun  8 22:11:20 2025]
rule minSBM:
    input: cache/trimmed_graph/fdr-1e-5/example.xml.gz
    output: cache/initialBlock/fdr-1e-5/example.dill
    log: logs/intialBlock/fdr-1e-5/example.log
    jobid: 8
    reason: Missing output files: cache/initialBlock/fdr-1e-5/example.dill; Input files updated by another job: cache/trimmed_graph/fdr-1e-5/example.xml.gz
    wildcards: fdr=1e-5, label=example
    threads: 8
    resources: tmpdir=/tmp

[Sun Jun  8 22:53:46 2025]
Finished job 8.
3 of 9 steps (33%) done
Select jobs to execute...

[Sun Jun  8 22:53:46 2025]
rule annealing:
    input: cache/trimmed_graph/fdr-1e-5/example.xml.gz, cache/initialBlock/fdr-1e-5/example.dill
    output: cache/annealedBlock/fdr-1e-5/example.dill
    log: logs/annealedBlock/fdr-1e-5/example.log
    jobid: 7
    reason: Missing output files: cache/annealedBlock/fdr-1e-5/example.dill; Input files updated by another job: cache/trimmed_graph/fdr-1e-5/example.xml.gz, cache/initialBlock/fdr-1e-5/example.dill
    wildcards: fdr=1e-5, label=example
    threads: 8
    resources: tmpdir=/tmp





[Mon Jun  9 10:36:54 2025]
Finished job 7.
4 of 9 steps (44%) done
Select jobs to execute...

[Mon Jun  9 10:36:54 2025]
rule equilibrate:
    input: cache/trimmed_graph/fdr-1e-5/example.xml.gz, cache/annealedBlock/fdr-1e-5/example.dill
    output: cache/equilibrate/fdr-1e-5/example.dill
    log: logs/equilibrate/fdr-1e-5/example.log
    jobid: 6
    reason: Missing output files: cache/equilibrate/fdr-1e-5/example.dill; Input files updated by another job: cache/trimmed_graph/fdr-1e-5/example.xml.gz, cache/annealedBlock/fdr-1e-5/example.dill
    wildcards: fdr=1e-5, label=example
    threads: 8
    resources: tmpdir=/tmp




