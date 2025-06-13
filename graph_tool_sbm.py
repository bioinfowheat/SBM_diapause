import graph_tool.all as gt
from graph_tool.all import *
import pandas as pd
import numpy as np
import dill
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm


### Load the graph
g = gt.load_graph("cache/trimmed_graph/fdr-1e-5/example.xml.gz")
state = gt.minimize_blockmodel_dl(g)


### Save the optimized SBM as a BlockState object
# Important to save the graph g AND the state from "cache/equilibrate/fdr-1e-5/example.dill" together
with open("cache/equilibrate/fdr-1e-5/example.dill", "wb") as f:
    dill.dump((g,state), f)


### Load the optimized SBM block state
with open("cache/equilibrate/fdr-1e-5/example.dill", "rb") as f:
    g, state = dill.load(f)


print(type(state)) #Should be BlockState object
blocks = state.get_blocks()  #returns a mapping of each vertex to its assigned block (at top level), based on the fitted (SBM) stored in the BlockState object state
print("block type =",type(blocks)) #Should be VertexPropertyMap


### Add block membership as vertex property (vp)
g.vp["block"] = blocks #adds a new vp called "block" to the graph g, and sets its values to the block IDs stored in blocks
print("list of vp =",list(g.vp.keys())) #double-check that "block" was added as a vp


### Load gene-block assignment
gene_blocks = pd.read_csv("cache/blockSummary/fdr-1e-5/example/gene_block.csv")  


### Load block summary
block_summary = pd.read_csv("cache/blockSummary/fdr-1e-5/example/block_summary.csv")


### Add gene names to the graph: assign it from gene_block.csv if order is correct:
gene_names = gene_blocks["Gene"].tolist()
g.vp["gene_name"] = g.new_vertex_property("string")
for v in g.vertices():
    g.vp["gene_name"][v] = gene_names[int(v)]


### Attach block size (from block_summary) as a vertex property
g.vp["block_size"] = g.new_vertex_property("int")
block_size_map = dict(zip(block_summary["Block"], block_summary["N_genes"]))


for v in g.vertices():
    blk = g.vp["block"][v]
    g.vp["block_size"][v] = block_size_map.get(blk, 0)


### Attach block average connectivity (from block_summary)
g.vp["block_connectivity"] = g.new_vertex_property("double")
connectivity_map = dict(zip(block_summary["Block"], block_summary["Average_total_degree"]))


for v in g.vertices():
    blk = g.vp["block"][v]
    g.vp["block_connectivity"][v] = connectivity_map.get(blk, 0)


### Edge widths based on weight available in our SBM object
print("Available edge properties:", list(g.ep.keys())) #to find out the names of the weights
eweight = g.ep["z_s"]


### Plot:
## Plot layout
pos = gt.sfdp_layout(g)


## Color vertices by block membership
blocks = g.vp["block"]
num_blocks = max(blocks.a) + 1
# Generate distinct colors for blocks using a colormap
cmap = plt.colormaps['gist_ncar'].resampled(num_blocks)
v_colors = [mcolors.to_rgba(cmap(b)) for b in blocks.a]
# Assign colors
g.vp["color"] = g.new_vertex_property("vector<float>")
for v in g.vertices():
    g.vp["color"][v] = v_colors[int(blocks[v])][:3]  # rgb only


## Set vertex sizes proportional to block_size (scale for visualization)
block_sizes = g.vp["block_size"].a
min_size = 5
max_size = 20
v_sizes_array = min_size + (block_sizes - block_sizes.min()) / (block_sizes.max() - block_sizes.min()) * (max_size - min_size)
# Convert vertex sizes numpy array to vertex property map
v_sizes = g.new_vertex_property("double")
for v in g.vertices():
    v_sizes[v] = v_sizes_array[int(v)]


## Set edge widths proportional to edge weight "z_s"
edge_weights = eweight.a
min_ew = 0.5
max_ew = 3
e_widths_array = min_ew + (edge_weights - edge_weights.min()) / (edge_weights.max() - edge_weights.min()) * (max_ew - min_ew)
# Convert edge widths numpy array to edge property map
e_widths = g.new_edge_property("double")
for e in g.edges():
    e_widths[e] = e_widths_array[g.edge_index[e]]


## Draw the graph
#Flat graph
gt.graph_draw(g,
              vertex_fill_color=g.vp["color"], #distinct color of vertices for each block
              vertex_size=v_sizes, #vertex size proportional to block size
              edge_pen_width=e_widths, #use edge weight for edge width
              output_size=(1000, 1000),
              output="graph_final.png")  # save to file


gt.graph_draw(g,
              vertex_fill_color=g.vp["color"], #distinct color of vertices for each block
              vertex_size=v_sizes, #vertex size proportional to block size
              edge_pen_width=e_widths, #use edge weight for edge width
              output_size=(1000, 1000),
              output="graph_final.svg")  # save to file