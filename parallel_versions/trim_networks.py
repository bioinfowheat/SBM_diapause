import os
# Set environment variables BEFORE any imports
os.environ["OMP_NUM_THREADS"] = "32"
os.environ["OPENBLAS_NUM_THREADS"] = "32"
os.environ["NUMEXPR_NUM_THREADS"] = "32"
os.environ["NUMBA_NUM_THREADS"] = "32"
os.environ["MKL_NUM_THREADS"] = "32"

import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception

from graph_tool.all import *
import pandas as pd
import numpy as np
import scipy as sp
from sklearn.covariance import LedoitWolf, OAS
import statsmodels.api as sm
from multiprocessing import Pool, cpu_count
from functools import partial

from multipy.fdr import lsu

def filterByEdge(g, corr, cutOff, keepOnlyMain):
    # Filtering edges
    corr = g.edge_properties[corr]
    sign = g.new_ep("bool", True)
    sign.a = np.array(np.abs(corr.a) > cutOff)

    tv = GraphView(g, efilt=sign)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

def parallel_fdr_correction(pvals_chunk, q_level):
    """Apply FDR correction to a chunk of p-values"""
    return lsu(pvals_chunk, q=q_level)

def filterByFDR_parallel(g, level, keepOnlyMain, n_cores=None):
    """
    Parallel FDR filtering with chunked processing
    """
    if n_cores is None:
        n_cores = min(32, cpu_count())
    
    logging.info(f"Using {n_cores} cores for FDR correction...")
    
    # Get p-values
    pvals = np.array(g.edge_properties["pvalue"].a)
    n_edges = len(pvals)
    
    if n_edges < 10000:  # For small graphs, don't parallelize
        logging.info("Small graph, using sequential FDR correction...")
        fdr_result = lsu(pvals, q=level)
    else:
        # Split p-values into chunks for parallel processing
        chunk_size = max(1000, n_edges // n_cores)
        chunks = [pvals[i:i + chunk_size] for i in range(0, n_edges, chunk_size)]
        
        # Process chunks in parallel
        with Pool(n_cores) as pool:
            fdr_func = partial(parallel_fdr_correction, q_level=level)
            chunk_results = pool.map(fdr_func, chunks)
        
        # Combine results
        fdr_result = np.concatenate(chunk_results)
    
    # Create edge filter
    fdr_ep = g.new_ep("bool", True)
    fdr_ep.a = fdr_result

    tv = GraphView(g, efilt=fdr_ep)

    # Keeping largest component
    if keepOnlyMain:
        logging.info("Finding largest connected component...")
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

def filterByFDR(g, level, keepOnlyMain):
    """Original sequential version for compatibility"""
    # Filtering edges
    pvals = np.array(g.edge_properties["pvalue"].a)

    fdr_ep = g.new_ep("bool", True)
    fdr_ep.a = lsu(pvals, q=level)

    tv = GraphView(g, efilt=fdr_ep)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

def parallel_arctanh_transform(spearman_chunk):
    """Apply arctanh transformation to correlation chunk"""
    return 2 * np.arctanh(spearman_chunk)

def add_z_weights_parallel(gi, n_cores=None):
    """
    Add z-transformed weights in parallel for large graphs
    """
    if n_cores is None:
        n_cores = min(32, cpu_count())
    
    spearman = gi.edge_properties["spearman"]
    spearman_array = spearman.a
    n_edges = len(spearman_array)
    
    if n_edges < 10000:  # For small graphs, use sequential
        logging.info("Small graph, using sequential z-transformation...")
        z_values = 2 * np.arctanh(spearman_array)
    else:
        logging.info(f"Using {n_cores} cores for z-transformation...")
        
        # Split into chunks
        chunk_size = max(1000, n_edges // n_cores)
        chunks = [spearman_array[i:i + chunk_size] for i in range(0, n_edges, chunk_size)]
        
        # Process in parallel
        with Pool(n_cores) as pool:
            chunk_results = pool.map(parallel_arctanh_transform, chunks)
        
        # Combine results
        z_values = np.concatenate(chunk_results)
    
    # Create new edge property
    gi.edge_properties["z_s"] = gi.new_edge_property("double", z_values)

def saveGeneTable(g, gene_expr, output=None, transpose = False):
    gene_list = []
    for i in g.vertex_properties['genes']:
        gene_list.append(i)
    expr_df = gene_expr[gene_list]
    if transpose is True:
        expr_df = expr_df.T
    if output is not None:
        expr_df.to_csv(output)
    return expr_df

if __name__ == '__main__':

    # Check if graph-tool has OpenMP support
    logging.info(f"Graph-tool OpenMP enabled: {openmp_enabled()}")
    logging.info(f"Available CPU cores: {cpu_count()}")
    
    logging.info("FDR level:" + str(snakemake.wildcards.fdr))
    fdr = float(snakemake.wildcards.fdr)

    logging.info("Reading full graph...")
    g = load_graph(snakemake.input[0])
    
    initial_edges = g.num_edges()
    initial_vertices = g.num_vertices()
    logging.info(f"Initial graph: {initial_vertices} vertices, {initial_edges} edges")

    logging.info("Trimming with parallel FDR correction...")
    
    # Use parallel version for large graphs, sequential for small ones
    if initial_edges > 50000:
        tv = filterByFDR_parallel(g, fdr, True, n_cores=32)
    else:
        tv = filterByFDR(g, fdr, True)
    
    logging.info("Creating pruned graph...")
    gi = Graph(tv, prune = True)

    N = len(gi.get_vertices())
    logging.info(str(N) + " genes")
    Et = (N * N - N)/2
    E = len(gi.get_edges())
    density = E/Et
    logging.info("Density: " + str(density))
    logging.info("Min correlation: " + str(min(gi.edge_properties["spearman"].a)))
    logging.info("Max correlation: " + str(max(gi.edge_properties["spearman"].a)))
    
    logging.info("Adding z-transformed weights...")
    add_z_weights_parallel(gi, n_cores=32)

    logging.info("Writing trimmed graph...")
    gi.save(snakemake.output[0])
    logging.info("Done!")