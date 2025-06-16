import graph_tool.all as gt
from graph_tool.all import *
import pandas as pd
import numpy as np
import dill
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import sys
import os


def process_graph_sbm_hierarchy(filename, fdr_threshold="1e-5"):
    """
    Process graph SBM hierarchy with configurable filename and FDR threshold.
    
    Args:
        filename: Base filename (without extension) for the graph files
        fdr_threshold: FDR threshold string (default: "1e-5")
    """
    
    # Construct file paths based on the input filename
    graph_path = f"cache/trimmed_graph/fdr-{fdr_threshold}/{filename}.xml.gz"
    state_save_path = f"cache/MCMC/blocks/fdr-{fdr_threshold}/{filename}.dill"
    gene_blocks_path = f"cache/blockSummary/fdr-{fdr_threshold}/{filename}/gene_block.csv"
    block_summary_path = f"cache/blockSummary/fdr-{fdr_threshold}/{filename}/block_summary.csv"
    
    # Output filenames
    output_png = f"graph_{filename}_hierarchy.png"
    output_svg = f"graph_{filename}_hierarchy.svg"
    output_png2 = f"graph_{filename}.png"
    output_svg2 = f"graph_{filename}.svg"
    
    print(f"Processing hierarchical graph: {filename}")
    print(f"FDR threshold: {fdr_threshold}")
    
    ### Load the graph
    print(f"Loading graph from: {graph_path}")
    g = gt.load_graph(graph_path)
    state = gt.minimize_nested_blockmodel_dl(g)

    ### Save the optimized nested SBM as a BlockState object
    # Important to save the graph g AND the state together
    print(f"Saving hierarchical state to: {state_save_path}")
    os.makedirs(os.path.dirname(state_save_path), exist_ok=True)
    with open(state_save_path, "wb") as f:
        dill.dump((g, state), f)

    ### Load the optimized nested SBM block state
    with open(state_save_path, "rb") as f:
        g, state = dill.load(f)

    print(type(state))  # Should be NestedBlockState object
    blocks = state.levels[0].get_blocks()  # returns a mapping of each vertex to its assigned block (at top level)
    print("block type =", type(blocks))  # Should be VertexPropertyMap

    ### Add block membership as vertex property (vp)
    g.vp["block"] = blocks
    print("list of vp =", list(g.vp.keys()))

    ### Load gene-block assignment
    print(f"Loading gene blocks from: {gene_blocks_path}")
    gene_blocks = pd.read_csv(gene_blocks_path)

    ### Load block summary
    print(f"Loading block summary from: {block_summary_path}")
    block_summary = pd.read_csv(block_summary_path)

    ### Add gene names to the graph
    gene_names = gene_blocks["Gene"].tolist()
    g.vp["gene_name"] = g.new_vertex_property("string")
    for v in g.vertices():
        g.vp["gene_name"][v] = gene_names[int(v)]

    ### Attach block size as a vertex property
    g.vp["block_size"] = g.new_vertex_property("int")
    block_size_map = dict(zip(block_summary["Block"], block_summary["N_genes"]))

    for v in g.vertices():
        blk = g.vp["block"][v]
        g.vp["block_size"][v] = block_size_map.get(blk, 0)

    ### Attach block average connectivity
    g.vp["block_connectivity"] = g.new_vertex_property("double")
    connectivity_map = dict(zip(block_summary["Block"], block_summary["Average_total_degree"]))

    for v in g.vertices():
        blk = g.vp["block"][v]
        g.vp["block_connectivity"][v] = connectivity_map.get(blk, 0)

    ### Edge widths based on weight available in our SBM object
    print("Available edge properties:", list(g.ep.keys()))
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

    ## Set vertex sizes proportional to block_size
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

    ## Draw the hierarchical graph
    print(f"Saving hierarchical PNG to: {output_png}")
    gt.draw_hierarchy(state,
                      output_size=(2000, 2000),
                      output=output_png)

    print(f"Saving hierarchical SVG to: {output_svg}")
    gt.draw_hierarchy(state,
                      output_size=(2000, 2000),
                      output=output_svg)

    ## Draw the graph
    print(f"Saving PNG to: {output_png2}")
    gt.graph_draw(g,
                  vertex_fill_color=g.vp["color"],
                  vertex_size=v_sizes,
                  edge_pen_width=e_widths,
                  output_size=(1000, 1000),
                  output=output_png)

    print(f"Saving SVG to: {output_svg2}")
    gt.graph_draw(g,
                  vertex_fill_color=g.vp["color"],
                  vertex_size=v_sizes,
                  edge_pen_width=e_widths,
                  output_size=(1000, 1000),
                  output=output_svg)
   
    print(f"Processing complete for {filename}")
    
    print(f"Hierarchical processing complete for {filename}")
    
    # Print hierarchy information
    print(f"Number of hierarchy levels: {len(state.levels)}")
    for i, level in enumerate(state.levels):
        num_blocks_at_level = level.get_nonempty_B()
        print(f"Level {i}: {num_blocks_at_level} blocks")


def main():
    """Main function to handle command line arguments or interactive input."""
    if len(sys.argv) >= 2:
        filename = sys.argv[1]
        fdr_threshold = sys.argv[2] if len(sys.argv) >= 3 else "1e-5"
    else:
        # Interactive input if no command line arguments
        filename = input("Enter the filename (without extension): ").strip()
        if not filename:
            print("Error: Filename cannot be empty")
            return
        
        fdr_threshold = input("Enter FDR threshold (default: 1e-5): ").strip()
        if not fdr_threshold:
            fdr_threshold = "1e-5"
    
    try:
        process_graph_sbm_hierarchy(filename, fdr_threshold)
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        print("Please check that all required files exist in the expected locations.")
    except Exception as e:
        print(f"Error processing hierarchical graph: {e}")


if __name__ == "__main__":
    main()