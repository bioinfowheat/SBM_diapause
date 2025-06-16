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

import dill
from multiprocessing import Pool, cpu_count

from graph_tool.all import *
import numpy as np

def run_parallel_annealing(state_min, beta_range, niter, n_chains=8):
    """
    Run parallel annealing with multiple independent chains
    """
    logging.info(f"Running parallel annealing with {n_chains} chains...")
    
    # Create multiple independent states
    states = [state_min.copy() for _ in range(n_chains)]
    
    # Different annealing schedules for diversity
    beta_ranges = []
    for i in range(n_chains):
        # Vary the annealing schedule slightly between chains
        beta_min = beta_range[0] + i * 0.1
        beta_max = beta_range[1] + i * 0.5
        beta_ranges.append((beta_min, beta_max))
    
    logging.info(f"Beta ranges for chains: {beta_ranges}")
    
    def anneal_single_chain(args):
        state, chain_id, beta_r, chain_niter = args
        
        try:
            # Perform annealing on this chain
            S1 = state.entropy()
            mcmc_anneal(state, beta_range=beta_r, niter=chain_niter,
                       mcmc_equilibrate_args=dict(force_niter=10), verbose=False)
            S2 = state.entropy()
            
            improvement = S2 - S1
            logging.info(f"Chain {chain_id}: improvement = {improvement}")
            
            return state, S1, S2, improvement
            
        except Exception as e:
            logging.error(f"Chain {chain_id} failed: {e}")
            return None, float('inf'), float('inf'), 0
    
    # Prepare arguments for parallel execution
    chain_niter = max(100, niter // 2)  # Each chain runs for half the iterations
    args = [(states[i], i, beta_ranges[i], chain_niter) for i in range(n_chains)]
    
    try:
        # Run chains in parallel
        with Pool(n_chains) as pool:
            results = pool.map(anneal_single_chain, args)
        
        # Find the best result
        best_state = None
        best_improvement = float('-inf')
        best_S1 = None
        best_S2 = None
        
        for state, S1, S2, improvement in results:
            if state is not None and improvement > best_improvement:
                best_state = state
                best_improvement = improvement
                best_S1 = S1
                best_S2 = S2
        
        if best_state is not None:
            logging.info(f"Best improvement from parallel annealing: {best_improvement}")
            return best_state, best_S1, best_S2, best_improvement
        else:
            logging.error("All parallel annealing chains failed!")
            return None, None, None, None
            
    except Exception as e:
        logging.error(f"Parallel annealing failed: {e}")
        return None, None, None, None

def run_multi_stage_annealing(state_min, beta_range, niter, n_stages=4):
    """
    Run multi-stage annealing with different temperature schedules
    """
    logging.info(f"Running multi-stage annealing with {n_stages} stages...")
    
    current_state = state_min.copy()
    total_improvement = 0
    S_initial = current_state.entropy()
    
    # Divide iterations across stages
    stage_niter = max(100, niter // n_stages)
    
    for stage in range(n_stages):
        logging.info(f"Annealing stage {stage + 1}/{n_stages}...")
        
        # Adjust beta range for each stage
        stage_factor = (stage + 1) / n_stages
        stage_beta_min = beta_range[0] * (1 + stage * 0.2)
        stage_beta_max = beta_range[1] * (1 + stage * 0.1)
        stage_beta_range = (stage_beta_min, stage_beta_max)
        
        logging.info(f"Stage {stage + 1} beta range: {stage_beta_range}")
        
        try:
            S1 = current_state.entropy()
            mcmc_anneal(current_state, beta_range=stage_beta_range, niter=stage_niter,
                       mcmc_equilibrate_args=dict(force_niter=10), verbose=False)
            S2 = current_state.entropy()
            
            stage_improvement = S2 - S1
            total_improvement += stage_improvement
            
            logging.info(f"Stage {stage + 1} improvement: {stage_improvement}")
            
        except Exception as e:
            logging.error(f"Stage {stage + 1} failed: {e}")
            break
    
    S_final = current_state.entropy()
    logging.info(f"Total multi-stage improvement: {total_improvement}")
    
    return current_state, S_initial, S_final, total_improvement

def run_adaptive_annealing(state_min, beta_range, niter):
    """
    Run adaptive annealing with dynamic parameter adjustment
    """
    logging.info("Running adaptive annealing...")
    
    current_state = state_min.copy()
    S1 = current_state.entropy()
    
    # Adaptive parameters
    current_beta_range = beta_range
    current_niter = niter
    improvement_history = []
    
    # Multiple annealing rounds with adaptation
    total_improvement = 0
    n_rounds = 3
    
    for round_num in range(n_rounds):
        logging.info(f"Adaptive annealing round {round_num + 1}/{n_rounds}")
        
        round_niter = max(200, current_niter // n_rounds)
        
        try:
            S_before = current_state.entropy()
            mcmc_anneal(current_state, beta_range=current_beta_range, niter=round_niter,
                       mcmc_equilibrate_args=dict(force_niter=10), verbose=False)
            S_after = current_state.entropy()
            
            round_improvement = S_after - S_before
            improvement_history.append(round_improvement)
            total_improvement += round_improvement
            
            logging.info(f"Round {round_num + 1} improvement: {round_improvement}")
            
            # Adapt parameters based on improvement
            if round_improvement > 0:
                # If improving, increase intensity
                current_beta_range = (current_beta_range[0] * 0.9, current_beta_range[1] * 1.1)
            else:
                # If not improving much, try different range
                current_beta_range = (current_beta_range[0] * 1.1, current_beta_range[1] * 0.9)
            
            logging.info(f"Adapted beta range for next round: {current_beta_range}")
            
        except Exception as e:
            logging.error(f"Adaptive annealing round {round_num + 1} failed: {e}")
            break
    
    S2 = current_state.entropy()
    logging.info(f"Total adaptive improvement: {total_improvement}")
    
    return current_state, S1, S2, total_improvement

def run_sequential_annealing(state_min, beta_range, niter):
    """
    Original sequential annealing for comparison/fallback
    """
    logging.info("Running sequential annealing...")
    
    S1 = state_min.entropy()
    mcmc_anneal(state_min, beta_range=beta_range, niter=niter,
                mcmc_equilibrate_args=dict(force_niter=10), verbose=False)
    S2 = state_min.entropy()
    
    improvement = S2 - S1
    logging.info(f"Sequential annealing improvement: {improvement}")
    
    return state_min, S1, S2, improvement

if __name__ == '__main__':

    # Check system capabilities
    logging.info(f"Graph-tool OpenMP enabled: {openmp_enabled()}")
    logging.info(f"Available CPU cores: {cpu_count()}")
    n_cores = min(8, cpu_count())

    logging.info("Reading trimmed graph...")
    g = load_graph(snakemake.input.graph)
    logging.info(f"Graph has {g.num_vertices()} vertices and {g.num_edges()} edges")
    
    logging.info("Loading blocks...")
    with open (snakemake.input.blocks, "rb") as fh:
        bs = dill.load(fh)
    bs[0] = g.own_property(bs[0])

    logging.info("Creating nested block model...")
    if snakemake.params.type == "layer": 
        state_min = NestedBlockState(g, bs=bs,
                                     state_args=dict(base_type=LayeredBlockState,
                                      ec=g.ep.layer, layers=True,
                                      recs=[g.ep.z_s], 
                                      rec_types=["real-normal"]))
    elif snakemake.params.type == "batch":
        state_min = NestedBlockState(g, bs=bs,
                                     state_args=dict(recs=[g.ep.z_s],
                                                     rec_types=["real-normal"]))
    else: 
        sys.exit('Parameter "type" is not batch or layer')

    # Annealing parameters
    beta_range = (1, 10)
    niter = 1000
    n_vertices = g.num_vertices()
    
    logging.info("Starting annealing...")
    
    # Choose annealing strategy based on graph size and available cores
    if n_vertices > 1000 and n_cores >= 4:
        # Strategy 1: Try parallel annealing for larger graphs
        logging.info("Attempting parallel annealing...")
        result = run_parallel_annealing(state_min, beta_range, niter, n_chains=n_cores)
        
        if result[0] is not None:
            final_state, S1, S2, improvement = result
            logging.info("Successfully used parallel annealing")
        else:
            # Fallback to multi-stage
            logging.info("Parallel annealing failed, trying multi-stage annealing...")
            final_state, S1, S2, improvement = run_multi_stage_annealing(state_min, beta_range, niter)
    
    elif n_vertices > 500:
        # Strategy 2: Multi-stage annealing for medium graphs
        logging.info("Using multi-stage annealing...")
        final_state, S1, S2, improvement = run_multi_stage_annealing(state_min, beta_range, niter, n_stages=4)
    
    elif n_vertices > 200:
        # Strategy 3: Adaptive annealing for smaller graphs
        logging.info("Using adaptive annealing...")
        final_state, S1, S2, improvement = run_adaptive_annealing(state_min, beta_range, niter)
    
    else:
        # Strategy 4: Sequential annealing for very small graphs
        logging.info("Using sequential annealing (optimal for small graphs)...")
        final_state, S1, S2, improvement = run_sequential_annealing(state_min, beta_range, niter)

    # If we somehow don't have a final state, use the original
    if final_state is None:
        logging.warning("All annealing strategies failed, using original state")
        final_state = state_min
        S1 = final_state.entropy()
        S2 = S1
        improvement = 0

    logging.info("Improvement from annealing: " + str(improvement))
    logging.info("Final entropy after annealing: " + str(final_state.entropy()))

    logging.info("Saving blockstate...")
    block_state = final_state.get_bs()
    with open(snakemake.output.blocks, 'wb') as fh:
        dill.dump(block_state, fh, recurse=True)
    logging.info("Done!")