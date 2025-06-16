import os
# Set environment variables BEFORE any imports
os.environ["OMP_NUM_THREADS"] = '8'
os.environ["OPENBLAS_NUM_THREADS"] = '8'
os.environ["NUMEXPR_NUM_THREADS"] = '8'
os.environ["NUMBA_NUM_THREADS"] = '8'
os.environ["MKL_NUM_THREADS"] = '8'

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
import statsmodels.api as sm
from multiprocessing import Pool, cpu_count
from functools import partial
import itertools

def compute_correlations_chunk(args):
    """
    Compute Spearman correlations for a chunk of gene pairs
    """
    gene_pairs, ge_scale, gene_cols = args
    
    results = []
    for i, j in gene_pairs:
        try:
            spearman_r = sp.stats.spearmanr(ge_scale.iloc[:,i], ge_scale.iloc[:,j])
            results.append((i, j, spearman_r[0], spearman_r[1]))
        except Exception as e:
            # Handle potential numerical issues
            results.append((i, j, 0.0, 1.0))
    
    return results

def parallel_correlation_computation(ge_scale, n_cores=None):
    """
    Compute all pairwise Spearman correlations in parallel
    """
    if n_cores is None:
        n_cores = min(8, cpu_count())
    
    n_genes = ge_scale.shape[1]
    
    # Generate all unique pairs (i, j) where i > j
    all_pairs = [(i, j) for i in range(n_genes) for j in range(i)]
    n_pairs = len(all_pairs)
    
    logging.info(f"Computing {n_pairs} correlations using {n_cores} cores...")
    
    if n_pairs < 1000:  # For small datasets, don't parallelize
        logging.info("Small dataset, using sequential correlation computation...")
        return compute_correlations_chunk((all_pairs, ge_scale, ge_scale.columns))[0]
    
    # Split pairs into chunks for parallel processing
    chunk_size = max(100, n_pairs // (n_cores * 4))  # More chunks than cores for better load balancing
    chunks = [all_pairs[i:i + chunk_size] for i in range(0, n_pairs, chunk_size)]
    
    # Prepare arguments for parallel processing
    args = [(chunk, ge_scale, ge_scale.columns) for chunk in chunks]
    
    # Process chunks in parallel
    with Pool(n_cores) as pool:
        chunk_results = pool.map(compute_correlations_chunk, args)
    
    # Flatten results
    all_results = []
    for chunk_result in chunk_results:
        all_results.extend(chunk_result)
    
    return all_results

def parallel_arctanh_transform(correlations, n_cores=None):
    """
    Apply arctanh transformation to correlations in parallel
    """
    if n_cores is None:
        n_cores = min(8, cpu_count())
    
    if len(correlations) < 10000:
        # Sequential for small datasets
        return [2 * np.arctanh(corr) for corr in correlations]
    
    # Split correlations into chunks
    chunk_size = max(1000, len(correlations) // n_cores)
    chunks = [correlations[i:i + chunk_size] for i in range(0, len(correlations), chunk_size)]
    
    def transform_chunk(corr_chunk):
        return [2 * np.arctanh(corr) for corr in corr_chunk]
    
    with Pool(n_cores) as pool:
        chunk_results = pool.map(transform_chunk, chunks)
    
    # Flatten results
    all_transformed = []
    for chunk_result in chunk_results:
        all_transformed.extend(chunk_result)
    
    return all_transformed

def process_layer_parallel(layer_data, n_cores=None):
    """
    Process a single layer with parallel correlation computation
    """
    k, ge, ge_scale, input_name, param_type = layer_data
    
    logging.info(f"Processing layer {k} - {input_name} with parallel computation...")
    
    n_genes = ge.shape[1]
    
    # Compute all correlations in parallel
    correlation_results = parallel_correlation_computation(ge_scale, n_cores)
    
    # Process results
    layer_edges = []
    unit_corr_layer = {}
    
    for i, j, spearman_corr, p_value in correlation_results:
        gene_i = ge.columns[i]
        gene_j = ge.columns[j]
        
        edge_data = {
            'gene_i': gene_i,
            'gene_j': gene_j,
            'spearman': spearman_corr,
            'pvalue': p_value,
            'layer': k if param_type == "layer" else None,
            'dataset': input_name if param_type == "layer" else None
        }
        
        layer_edges.append(edge_data)
        
        # Check for unit correlations
        if abs(spearman_corr) > 0.999:
            logging.info(f"Warning: Spearman correlation between {gene_i} and {gene_j} is near unity. Flagging for removal.")
            if gene_i not in unit_corr_layer and gene_j not in unit_corr_layer:
                unit_corr_layer[gene_i] = True
    
    logging.info(f"Finished layer {k} with {len(layer_edges)} edges.")
    return layer_edges, unit_corr_layer

if __name__ == '__main__':
    
    # Check system capabilities
    logging.info(f"Available CPU cores: {cpu_count()}")
    n_cores = min(8, cpu_count())
    
    # Load data based on type
    if snakemake.params.type == "layer": 
        logging.info("Building layered network...")
        data_dir = snakemake.input[0]
        input_list = []
        for file in os.listdir(data_dir):
            if file.endswith(".csv"):
                input_list.append(file)
        input_names = list(map(lambda p: p[:p.rfind('.')], input_list))
        gene_expr = []
        for file in input_list:
            gene_expr_raw = pd.read_table(os.path.join(data_dir, file))
            gene_expr.append(gene_expr_raw.T)
    elif snakemake.params.type == "batch":
        logging.info("Building network...")
        gene_expr_raw = pd.read_table(snakemake.input[0])
        gene_expr = [gene_expr_raw.T]
        input_names = ["batch"]
    else: 
        sys.exit('Parameter "type" is not batch or layer')
    
    logging.info("Read expression files...")

    n_layers = len(gene_expr)
    logging.info(f"Processing {n_layers} layers...")

    # Parallel data centering and scaling
    logging.info("Centering and scaling data matrices...")
    X_centered = []
    
    def center_and_scale(ge):
        return (ge - ge.mean()) / np.sqrt(ge.var())
    
    if n_layers > 1 and n_layers <= n_cores:
        # Parallelize across layers if we have multiple layers but not too many
        with Pool(min(n_layers, n_cores)) as pool:
            X_centered = pool.map(center_and_scale, gene_expr)
    else:
        # Sequential processing for single layer or too many layers
        for ge in gene_expr:
            X_centered.append(center_and_scale(ge))

    # Collect all gene names
    gene_names = set()
    number_of_genes = ""
    for ge in gene_expr:
        gene_names = gene_names.union(set(ge.columns))
        number_of_genes = number_of_genes + ", " + str(len(ge.columns)) + " genes"
    gene_names = list(gene_names)
    logging.info("Concatenating gene names" + number_of_genes)

    n_genes = len(gene_names)
    logging.info("Total number of genes: " + str(n_genes))

    # Estimate computational complexity
    total_pairs = sum(len(ge.columns) * (len(ge.columns) - 1) // 2 for ge in gene_expr)
    logging.info(f"Total correlation pairs to compute: {total_pairs}")

    # Create graph
    logging.info("Creating graph...")
    g = Graph(directed=False)
    g.add_vertex(n=n_genes)
    genes = g.new_vertex_property("string", np.array(gene_names, dtype="str"))
    g.vertex_properties["genes"] = genes

    # Initialize edge properties
    spearman = g.new_ep("double", 0)
    pval = g.new_ep("double", 0)
    layer = g.new_ep("int", 0)
    dataset = g.new_ep("string", "")

    # Process layers with parallel correlation computation
    all_edges = []
    unit_corr = {}
    
    # Determine processing strategy
    if n_layers == 1 or total_pairs > 100000:
        # For single layer or very large datasets, process sequentially but with parallel correlations
        for k in range(n_layers):
            layer_data = (k, gene_expr[k], X_centered[k], input_names[k], snakemake.params.type)
            edges, unit_corr_layer = process_layer_parallel(layer_data, n_cores)
            all_edges.extend(edges)
            unit_corr.update(unit_corr_layer)
    else:
        # For multiple smaller layers, process layers in parallel
        logging.info(f"Processing {n_layers} layers in parallel...")
        layer_args = [(k, gene_expr[k], X_centered[k], input_names[k], snakemake.params.type) 
                      for k in range(n_layers)]
        
        with Pool(min(n_layers, n_cores)) as pool:
            layer_results = pool.starmap(process_layer_parallel, 
                                       [(args, max(1, n_cores // n_layers)) for args in layer_args])
        
        # Combine results
        for edges, unit_corr_layer in layer_results:
            all_edges.extend(edges)
            unit_corr.update(unit_corr_layer)

    # Add edges to graph
    logging.info(f"Adding {len(all_edges)} edges to graph...")
    
    # Create gene name to vertex mapping for faster lookup
    gene_to_vertex = {}
    for v in g.vertices():
        gene_to_vertex[g.vp.genes[v]] = v
    
    for edge_data in all_edges:
        try:
            v_i = gene_to_vertex[edge_data['gene_i']]
            v_j = gene_to_vertex[edge_data['gene_j']]
            
            e = g.add_edge(v_i, v_j)
            spearman[e] = edge_data['spearman']
            pval[e] = edge_data['pvalue']
            
            if snakemake.params.type == "layer":
                layer[e] = edge_data['layer']
                dataset[e] = edge_data['dataset']
                
        except KeyError as ke:
            logging.warning(f"Gene not found in vertex mapping: {ke}")
            continue

    logging.info("Setting edge properties...")
    g.edge_properties["pvalue"] = pval
    g.edge_properties["spearman"] = spearman    
    
    if snakemake.params.type == "layer":
        logging.info("Setting layer properties...")
        g.edge_properties["layer"] = layer
        g.edge_properties["dataset"] = dataset
    
    # Remove vertices with unit correlations
    vertices_to_remove = []
    for gene in unit_corr.keys():
        logging.info("Removing vertex " + str(gene) + "!")
        try:
            v = gene_to_vertex[gene]
            vertices_to_remove.append(v)
        except KeyError:
            logging.warning(f"Gene {gene} not found for removal")
    
    if vertices_to_remove:
        g.remove_vertex(vertices_to_remove)
    
    # Parallel z-score transformation
    logging.info("Setting edge weights with parallel transformation...")
    correlations = g.ep.spearman.a
    
    if len(correlations) > 10000:
        logging.info("Using parallel arctanh transformation...")
        z_values = parallel_arctanh_transform(correlations, n_cores)
        g.edge_properties["z_s"] = g.new_edge_property("double", z_values)
    else:
        logging.info("Using sequential arctanh transformation...")
        g.edge_properties["z_s"] = g.new_edge_property("double", (2*np.arctanh(correlations)))
    
    logging.info("Writing graph...")
    g.save(snakemake.output[0])
    logging.info("Done!")