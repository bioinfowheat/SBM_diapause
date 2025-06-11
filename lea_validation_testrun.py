import graph_tool.all as gt
import pandas as pd
import numpy as np
import dill

### files needed are the following from an example.tsv run
### cache/blockSummary/fdr-1e-5/example/gene_block.csv
### cache/blockSummary/fdr-1e-5/example/block_summary.csv
### cache/trimmed_graph/fdr-1e-5/example.xml.gz
### cache/equilibrate/fdr-1e-5/example.dill

### Load the graph
g = gt.load_graph("cache/trimmed_graph/fdr-1e-5/example.xml.gz")

### Load the optimized SBM block state
with open("cache/equilibrate/fdr-1e-5/example.dill", "rb") as f:
    state = dill.load(f)

blocks = state.get_blocks()

### Add block membership as vertex property
g.vp["block"] = blocks.a.copy()  

### Load gene-block assignment
gene_blocks = pd.read_csv("cache/blockSummary/fdr-1e-5/example/gene_block.csv")  

### Load block summary 
block_summary = pd.read_csv("cache/blockSummary/fdr-1e-5/example/block_summary.csv") 

### Load GO annotation per gene (optional)
go_df = pd.read_csv("gene_go_annotations.csv") 

### Add gene names to the graph: assign it from gene_block.csv if order is correct:
gene_names = gene_blocks["gene"].tolist()
g.vp["gene_name"] = g.new_vertex_property("string")
for v in g.vertices():
    g.vp["gene_name"][v] = gene_names[int(v)]

### Attach GO annotations to vertices
g.vp["GO_terms"] = g.new_vertex_property("string")
go_dict = go_df.set_index("gene")["GO_terms"].to_dict()

for v in g.vertices():
    gene = g.vp["gene_name"][v]
    g.vp["GO_terms"][v] = go_dict.get(gene, "")

### Attach block size (from block_summary) as a vertex property
g.vp["block_size"] = g.new_vertex_property("int")
block_size_map = dict(zip(block_summary["block"], block_summary["size"]))

for v in g.vertices():
    blk = g.vp["block"][v]
    g.vp["block_size"][v] = block_size_map.get(blk, 0)

### Attach block average connectivity (from block_summary)
g.vp["block_connectivity"] = g.new_vertex_property("double")
connectivity_map = dict(zip(block_summary["block"], block_summary["avg_degree"]))

for v in g.vertices():
    blk = g.vp["block"][v]
    g.vp["block_connectivity"][v] = connectivity_map.get(blk, 0)
