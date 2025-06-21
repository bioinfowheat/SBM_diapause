
# getting output files for local plotting.
rawdatafile=normalized_counts
pval=fdr-1e-4
zip $rawdatafile.$pval.zip \
cache/blockSummary/$pval/$rawdatafile/gene_block.csv \
cache/blockSummary/$pval/$rawdatafile/block_summary.csv \
cache/trimmed_graph/$pval/$rawdatafile.xml.gz \
cache/MCMC/blocks/$pval/$rawdatafile.dill


# plotting when in the snakemake folder of a given run, 
# looking at results from a file called normalized_counts, where it was run with the config.yaml set to 1e-4 pvalue
python graph_tool_sbm_hiearchy.wheat_v1.py normalized_counts 1e-4
python graph_tool_sbm.wheat_v1.py normalized_counts 1e-4                                                                                  

# from server to a local computer, where I am in the local computer and running this to conenct back to the server to pull from
scp chrwhe@duke.zoologi.su.se:/mnt/griffin/chrwhe/SBM_2/test_run4/SBM-tools/snakemake/\*.zip .
scp chrwhe@duke.zoologi.su.se:/mnt/griffin/chrwhe/SBM_2/test_run4/SBM-tools/snakemake/\*.py .
scp chrwhe@duke.zoologi.su.se:/mnt/griffin/chrwhe/SBM_2/test_run4/SBM-tools/snakemake/\*.png .


# transfer
sudo cp normalized_counts_3.fdr-1e-2.zip /mnt/│griffin/leafar                 
sudo chown leafar:leafar /mnt/griffin/leafar/normalized_counts_3.fdr-1e-2.zip                                                                                                     

## 6a
