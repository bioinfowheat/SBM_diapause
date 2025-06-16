import os
# Set environment variables BEFORE any imports
os.environ["OMP_NUM_THREADS"] = '16'
os.environ["OPENBLAS_NUM_THREADS"] = '16'
os.environ["NUMEXPR_NUM_THREADS"] = '16'
os.environ["NUMBA_NUM_THREADS"] = '16'
os.environ["MKL_NUM_THREADS"] = '16'

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

import dill
from multiprocessing import Pool, cpu_count
from functools import partial
import concurrent.futures

from graph_tool.all import *
import numpy as np
import pandas as pd

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(file_path):
        os.makedirs(file_path)

def get_group(x, state):
    levels = state.get_levels()
    n_levels = len(levels)
    r = np.zeros(n_levels)
    r[0] = levels[0].get_blocks()[x]
    for i in range(1, n_levels):
        r[i] = levels[i].get_blocks()[r[i-1]]
    r = r.astype(int)
    return r

def parallel_get_group_chunk(args):
    """
    Process a chunk of vertices to get their group assignments
    """
    vertex_chunk, state = args
    
    results = []
    for v in vertex_chunk:
        group = get_group(v, state)
        results.append((v, group))
    
    return results

def create_nestedBlock_df_parallel(g, corr, state, n_cores=None):
    """
    Create nested block dataframe with parallel processing
    """
    if n_cores is None:
        n_cores = min(16, cpu_count())
    
    genes = g.vertex_properties["genes"]
    n_vertices = g.num_vertices()
    
    logging.info(f"Creating nested block DataFrame for {n_vertices} vertices using {n_cores} cores...")
    
    # Initialize DataFrame
    nested_block_df = pd.DataFrame(columns=('Gene', "Degree", "E_corr", 'B1', "B2", "B3", "B4", "B5", "B6", "B7"))
    
    vertices = list(g.vertex_index)
    
    if n_vertices < 1000:
        # Sequential processing for small graphs
        logging.info("Small graph, using sequential processing...")
        for v in vertices:
            line = [genes[v]]
            line.append(g.get_total_degrees([v])[0])
            line.append(np.mean(np.abs(g.get_all_edges(v, [corr])[:,2])))
            [line.append(i) for i in get_group(v, state)]
            nested_block_df.loc[v] = line
    else:
        # Parallel processing for large graphs
        
        # Step 1: Parallel computation of group assignments
        chunk_size = max(50, n_vertices // (n_cores * 2))
        vertex_chunks = [vertices[i:i + chunk_size] for i in range(0, n_vertices, chunk_size)]
        
        logging.info(f"Processing group assignments in {len(vertex_chunks)} chunks...")
        
        args = [(chunk, state) for chunk in vertex_chunks]
        
        with Pool(n_cores) as pool:
            chunk_results = pool.map(parallel_get_group_chunk, args)
        
        # Flatten results
        group_assignments = {}
        for chunk_result in chunk_results:
            for v, group in chunk_result:
                group_assignments[v] = group
        
        # Step 2: Parallel computation of other vertex properties
        def compute_vertex_properties_chunk(vertex_chunk):
            chunk_results = []
            for v in vertex_chunk:
                try:
                    gene_name = genes[v]
                    degree = g.get_total_degrees([v])[0]
                    
                    # Get edges and correlations for this vertex
                    edges_data = g.get_all_edges(v, [corr])
                    if len(edges_data) > 0:
                        e_corr = np.mean(np.abs(edges_data[:, 2]))
                    else:
                        e_corr = 0.0
                    
                    chunk_results.append((v, gene_name, degree, e_corr))
                except Exception as e:
                    logging.warning(f"Error processing vertex {v}: {e}")
                    chunk_results.append((v, f"vertex_{v}", 0, 0.0))
            
            return chunk_results
        
        logging.info("Computing vertex properties in parallel...")
        
        with Pool(n_cores) as pool:
            property_results = pool.map(compute_vertex_properties_chunk, vertex_chunks)
        
        # Combine all results into DataFrame
        logging.info("Assembling final DataFrame...")
        
        for chunk_result in property_results:
            for v, gene_name, degree, e_corr in chunk_result:
                line = [gene_name, degree, e_corr]
                [line.append(i) for i in group_assignments[v]]
                nested_block_df.loc[v] = line
    
    return nested_block_df

def parallel_process_block_level(args):
    """
    Process a single block level in parallel
    """
    level_data, block_dir, ers_data = args
    i, bl, block_df, n_levels, state = level_data
    
    level_results = []
    
    logging.info(f"Processing level {i+1} with {len(bl)} blocks...")
    
    for b in bl:
        try:
            line = [i+1, b]
            
            df = block_df[block_df['B' + str(i+1)] == b]
            genes = df["Gene"]
            file_name = "/" + '-'.join([str(num) for num in list(df.filter(like='B', axis=1).iloc[0, range(i, n_levels)])]) + ".csv"
            
            line.append(file_name)
            N_genes = genes.shape[0]
            line.append(N_genes)
            
            # Ensure directory exists and save genes
            ensure_dir(block_dir + "/Level_" + str(i+1))
            genes.to_csv(block_dir + "/Level_" + str(i+1) + file_name, header=False, index=False)
            
            # Calculate weighted metrics
            ers = ers_data[i]
            B = len(bl)
            E = ers.sum()
            
            if E > 0:
                q_r = (B/E) * (ers[b,b] - (ers[b,:].sum()**2/E))
                if np.abs(q_r) > 1:
                    logging.warning(f"q_r is larger than one for level {i+1}, block {b}: {q_r}")
                    q_r = np.clip(q_r, -1, 1)
            else:
                q_r = 0.0
            
            line.append(ers[b,b])
            line.append(ers[b,b] / N_genes if N_genes > 0 else 0)
            line.append(ers[b,:].sum())
            line.append(ers[b,:].sum() / N_genes if N_genes > 0 else 0)
            line.append(q_r)
            
            level_results.append(line)
            
        except Exception as e:
            logging.error(f"Error processing block {b} at level {i+1}: {e}")
            continue
    
    return level_results

def parallel_compute_ers_matrices(state, n_cores=None):
    """
    Compute ERS matrices for all levels in parallel
    """
    if n_cores is None:
        n_cores = min(16, cpu_count())
    
    levels = state.get_levels()
    n_levels = len(levels)
    
    logging.info(f"Computing ERS matrices for {n_levels} levels...")
    
    def compute_ers_for_level(level_state):
        try:
            return adjacency(level_state.bg, weight=level_state.mrs)
        except Exception as e:
            logging.error(f"Error computing ERS matrix: {e}")
            return None
    
    if n_levels <= n_cores:
        # Parallelize across levels
        with Pool(min(n_levels, n_cores)) as pool:
            ers_matrices = pool.map(compute_ers_for_level, levels)
    else:
        # Sequential computation if too many levels
        ers_matrices = [compute_ers_for_level(level) for level in levels]
    
    return ers_matrices

if __name__ == '__main__':

    # Check system capabilities
    logging.info(f"Available CPU cores: {cpu_count()}")
    n_cores = min(16, cpu_count())

    logging.info("Reading trimmed graph...")
    g = load_graph(snakemake.input.graph)
    logging.info(f"Graph has {g.num_vertices()} vertices and {g.num_edges()} edges")

    logging.info("Arctan transform on correlations...")
    corr = g.edge_properties["spearman"]
    
    # Parallel arctanh transformation for large graphs
    if g.num_edges() > 10000:
        logging.info("Using parallel arctanh transformation...")
        
        def parallel_arctanh(corr_chunk):
            return 2 * np.arctanh(corr_chunk)
        
        corr_array = corr.a
        chunk_size = max(1000, len(corr_array) // n_cores)
        chunks = [corr_array[i:i + chunk_size] for i in range(0, len(corr_array), chunk_size)]
        
        with Pool(n_cores) as pool:
            transformed_chunks = pool.map(parallel_arctanh, chunks)
        
        z_values = np.concatenate(transformed_chunks)
        g.ep.z_s = g.new_edge_property("double", z_values)
    else:
        logging.info("Using sequential arctanh transformation...")
        g.ep.z_s = g.new_edge_property("double", (2*np.arctanh(corr.a)))
    
    logging.info("Loading blocks...")
    with open (snakemake.input.blocks, "rb") as fh:
        bs = dill.load(fh)[0:7]
    bs[0] = g.own_property(bs[0])

    logging.info("Creating nested block model...")
    if snakemake.params.type == "layer": 
        state = NestedBlockState(g, bs=bs,
                                     state_args=dict(base_type=LayeredBlockState,
                                      ec=g.ep.layer, layers=True,
                                      recs=[g.ep.z_s], 
                                      rec_types=["real-normal"]))
    elif snakemake.params.type == "batch":
        state = NestedBlockState(g, bs=bs,
                                     state_args=dict(recs=[g.ep.z_s],
                                                     rec_types=["real-normal"]))
    else: 
        sys.exit('Parameter "type" is not batch or layer')

    logging.info("Creating block DataFrame with parallel processing...")
    block_df = create_nestedBlock_df_parallel(g, corr, state, n_cores)

    logging.info("Calculating block sizes...")
    blocks = [list(set(block_df[b])) for b in block_df.filter(like='B', axis=1)]
    block_sizes = [len(b) for b in blocks]
    block_sizes = -np.sort(-np.array(list(set(block_sizes))))
    block_sizes = [x for x in block_sizes if x >= 2]

    # Save background genes
    block_df["Gene"].to_csv(snakemake.params.blockDir + "/background.csv", header=False, index=False)

    logging.info("Pre-computing ERS matrices...")
    ers_matrices = parallel_compute_ers_matrices(state, n_cores)

    logging.info("Creating gene lists with parallel processing...")
    n_levels = len(block_sizes)
    output_df = pd.DataFrame(columns=('Nested_Level', 'Block', 'File', 'N_genes', 'Internal_degree', 'Average_internal_degree', 'Total_degree', 'Average_total_degree', 'Assortativity'))
    
    # Determine processing strategy
    total_blocks = sum(len(bl) for bl in blocks)
    
    if total_blocks > 100 and n_cores >= 4:
        # Parallel processing for many blocks
        logging.info(f"Processing {total_blocks} blocks across {n_levels} levels in parallel...")
        
        # Prepare arguments for parallel processing
        level_args = []
        for i in range(n_levels):
            level_data = (i, blocks[i], block_df, n_levels, state)
            level_args.append((level_data, snakemake.params.blockDir, ers_matrices))
        
        # Process levels in parallel
        with Pool(min(n_levels, n_cores)) as pool:
            all_level_results = pool.map(parallel_process_block_level, level_args)
        
        # Combine results
        l = 0
        for level_results in all_level_results:
            for line in level_results:
                output_df.loc[l] = line
                l += 1
    
    else:
        # Sequential processing for small datasets
        logging.info("Using sequential processing for block analysis...")
        l = 0
        for i in range(n_levels):
            logging.info("At level: " + str(i+1))
            bl = blocks[i]
            for b in bl:
                line = [i+1, b]
                
                df = block_df[block_df['B' + str(i+1)] == b]
                genes = df["Gene"]
                file_name = "/" + '-'.join([str(num) for num in list(df.filter(like='B', axis=1).iloc[0, range(i, n_levels)])]) + ".csv"
                
                line.append(file_name)
                N_genes = genes.shape[0]
                line.append(N_genes)
                
                ensure_dir(snakemake.params.blockDir + "/Level_" + str(i+1))
                genes.to_csv(snakemake.params.blockDir + "/Level_" + str(i+1) + file_name, header=False, index=False)
                
                # Weighted calculations
                ers = ers_matrices[i]
                B = len(bl)
                E = ers.sum()
                
                if E > 0:
                    q_r = (B/E) * (ers[b,b] - (ers[b,:].sum()**2/E))
                    if np.abs(q_r) > 1:
                        logging.warning(f"q_r is larger than one: {q_r}")
                        q_r = np.clip(q_r, -1, 1)
                else:
                    q_r = 0.0
                
                line.append(ers[b,b])
                line.append(ers[b,b] / N_genes if N_genes > 0 else 0)
                line.append(ers[b,:].sum())
                line.append(ers[b,:].sum() / N_genes if N_genes > 0 else 0)
                line.append(q_r)
                
                output_df.loc[l] = line
                l += 1

    logging.info("Outputting block summary...")
    output_df.to_csv(snakemake.output.blockSummary, index=False)

    logging.info("Outputting block DataFrame...")
    block_df.to_csv(snakemake.output.blockDF, index=False)

    logging.info("Done!")