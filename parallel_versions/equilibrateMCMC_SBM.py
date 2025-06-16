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
import time

from graph_tool.all import *
import numpy as np

def run_parallel_equilibration(state_min, wait_time, n_chains=8):
    """
    Run parallel MCMC equilibration with multiple independent chains
    """
    logging.info(f"Running parallel equilibration with {n_chains} chains...")
    
    # Create multiple independent states
    states = [state_min.copy() for _ in range(n_chains)]
    
    # Different wait times for diversity (some chains may converge faster)
    wait_times = []
    for i in range(n_chains):
        # Vary wait times: some shorter, some longer
        factor = 0.7 + (i * 0.1)  # Range from 0.7x to 1.4x
        chain_wait = max(100, int(wait_time * factor))
        wait_times.append(chain_wait)
    
    logging.info(f"Wait times for chains: {wait_times}")
    
    def equilibrate_single_chain(args):
        state, chain_id, chain_wait = args
        
        try:
            start_time = time.time()
            S1 = state.entropy()
            
            # Perform equilibration on this chain
            mcmc_equilibrate(state, wait=chain_wait, mcmc_args=dict(niter=10), verbose=False)
            
            S2 = state.entropy()
            end_time = time.time()
            
            improvement = S2 - S1
            runtime = end_time - start_time
            
            logging.info(f"Chain {chain_id}: improvement = {improvement:.6f}, runtime = {runtime:.2f}s")
            
            return state, S1, S2, improvement, runtime, True
            
        except Exception as e:
            logging.error(f"Chain {chain_id} failed: {e}")
            return None, float('inf'), float('inf'), 0, 0, False
    
    # Prepare arguments for parallel execution
    args = [(states[i], i, wait_times[i]) for i in range(n_chains)]
    
    try:
        # Run chains in parallel
        start_total = time.time()
        with Pool(n_chains) as pool:
            results = pool.map(equilibrate_single_chain, args)
        end_total = time.time()
        
        # Analyze results
        successful_results = [(state, S1, S2, improvement, runtime) 
                            for state, S1, S2, improvement, runtime, success in results if success]
        
        if not successful_results:
            logging.error("All parallel equilibration chains failed!")
            return None, None, None, None
        
        # Find the best result (lowest final entropy)
        best_state = None
        best_S2 = float('inf')
        best_S1 = None
        best_improvement = None
        
        for state, S1, S2, improvement, runtime in successful_results:
            if S2 < best_S2:
                best_state = state
                best_S2 = S2
                best_S1 = S1
                best_improvement = improvement
        
        logging.info(f"Parallel equilibration completed in {end_total - start_total:.2f}s")
        logging.info(f"Best improvement from parallel equilibration: {best_improvement:.6f}")
        logging.info(f"Successful chains: {len(successful_results)}/{n_chains}")
        
        return best_state, best_S1, best_S2, best_improvement
        
    except Exception as e:
        logging.error(f"Parallel equilibration failed: {e}")
        return None, None, None, None

def run_staged_equilibration(state_min, wait_time, n_stages=3):
    """
    Run multi-stage equilibration with progressively refined parameters
    """
    logging.info(f"Running staged equilibration with {n_stages} stages...")
    
    current_state = state_min.copy()
    total_improvement = 0
    S_initial = current_state.entropy()
    
    # Divide wait time across stages with increasing precision
    stage_waits = []
    for stage in range(n_stages):
        # Earlier stages use shorter waits, later stages longer
        stage_factor = 0.3 + (stage * 0.35)  # 0.3, 0.65, 1.0
        stage_wait = max(50, int(wait_time * stage_factor))
        stage_waits.append(stage_wait)
    
    logging.info(f"Stage wait times: {stage_waits}")
    
    for stage in range(n_stages):
        logging.info(f"Equilibration stage {stage + 1}/{n_stages} (wait={stage_waits[stage]})...")
        
        try:
            start_time = time.time()
            S1 = current_state.entropy()
            
            # Different niter for different stages
            stage_niter = 5 + (stage * 5)  # 5, 10, 15
            
            mcmc_equilibrate(current_state, wait=stage_waits[stage], 
                           mcmc_args=dict(niter=stage_niter), verbose=False)
            
            S2 = current_state.entropy()
            end_time = time.time()
            
            stage_improvement = S2 - S1
            total_improvement += stage_improvement
            
            logging.info(f"Stage {stage + 1}: improvement = {stage_improvement:.6f}, "
                        f"runtime = {end_time - start_time:.2f}s")
            
        except Exception as e:
            logging.error(f"Equilibration stage {stage + 1} failed: {e}")
            break
    
    S_final = current_state.entropy()
    logging.info(f"Total staged improvement: {total_improvement:.6f}")
    
    return current_state, S_initial, S_final, total_improvement

def run_adaptive_equilibration(state_min, wait_time):
    """
    Run adaptive equilibration with dynamic convergence detection
    """
    logging.info("Running adaptive equilibration...")
    
    current_state = state_min.copy()
    S1 = current_state.entropy()
    
    # Adaptive parameters
    current_wait = max(100, wait_time // 3)
    max_attempts = 5
    convergence_threshold = 1e-6
    
    total_improvement = 0
    
    for attempt in range(max_attempts):
        logging.info(f"Adaptive equilibration attempt {attempt + 1}/{max_attempts} (wait={current_wait})...")
        
        try:
            start_time = time.time()
            S_before = current_state.entropy()
            
            mcmc_equilibrate(current_state, wait=current_wait, 
                           mcmc_args=dict(niter=10), verbose=False)
            
            S_after = current_state.entropy()
            end_time = time.time()
            
            attempt_improvement = S_after - S_before
            total_improvement += attempt_improvement
            
            logging.info(f"Attempt {attempt + 1}: improvement = {attempt_improvement:.6f}, "
                        f"runtime = {end_time - start_time:.2f}s")
            
            # Check for convergence
            if abs(attempt_improvement) < convergence_threshold:
                logging.info(f"Convergence detected after attempt {attempt + 1}")
                break
            
            # Adapt wait time based on improvement
            if attempt_improvement > 0:
                # If improving, try longer equilibration
                current_wait = min(wait_time, int(current_wait * 1.3))
            else:
                # If not improving, try shorter equilibration or stop
                if attempt_improvement < -convergence_threshold:
                    logging.info("Deterioration detected, stopping adaptive equilibration")
                    break
                current_wait = max(50, int(current_wait * 0.8))
            
            logging.info(f"Adapted wait time for next attempt: {current_wait}")
            
        except Exception as e:
            logging.error(f"Adaptive equilibration attempt {attempt + 1} failed: {e}")
            break
    
    S2 = current_state.entropy()
    logging.info(f"Total adaptive improvement: {total_improvement:.6f}")
    
    return current_state, S1, S2, total_improvement

def run_sequential_equilibration(state_min, wait_time):
    """
    Original sequential equilibration for comparison/fallback
    """
    logging.info("Running sequential equilibration...")
    
    start_time = time.time()
    S1 = state_min.entropy()
    
    mcmc_equilibrate(state_min, wait=wait_time, mcmc_args=dict(niter=10), verbose=False)
    
    S2 = state_min.entropy()
    end_time = time.time()
    
    improvement = S2 - S1
    logging.info(f"Sequential equilibration: improvement = {improvement:.6f}, "
                f"runtime = {end_time - start_time:.2f}s")
    
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

    # Equilibration parameters
    wait_time = snakemake.params.wait
    n_vertices = g.num_vertices()
    
    logging.info(f"Starting MCMC equilibration with wait time: {wait_time}...")
    
    # Choose equilibration strategy based on graph size, wait time, and available cores
    if n_vertices > 1000 and wait_time > 500 and n_cores >= 4:
        # Strategy 1: Parallel equilibration for large graphs with long wait times
        logging.info("Attempting parallel equilibration...")
        result = run_parallel_equilibration(state_min, wait_time, n_chains=n_cores)
        
        if result[0] is not None:
            final_state, S1, S2, improvement = result
            logging.info("Successfully used parallel equilibration")
        else:
            # Fallback to staged equilibration
            logging.info("Parallel equilibration failed, trying staged equilibration...")
            final_state, S1, S2, improvement = run_staged_equilibration(state_min, wait_time)
    
    elif n_vertices > 500 or wait_time > 1000:
        # Strategy 2: Staged equilibration for medium graphs or long wait times
        logging.info("Using staged equilibration...")
        final_state, S1, S2, improvement = run_staged_equilibration(state_min, wait_time, n_stages=3)
    
    elif n_vertices > 200 and wait_time > 200:
        # Strategy 3: Adaptive equilibration for smaller graphs
        logging.info("Using adaptive equilibration...")
        final_state, S1, S2, improvement = run_adaptive_equilibration(state_min, wait_time)
    
    else:
        # Strategy 4: Sequential equilibration for small graphs or short wait times
        logging.info("Using sequential equilibration (optimal for small graphs/short waits)...")
        final_state, S1, S2, improvement = run_sequential_equilibration(state_min, wait_time)

    # If we somehow don't have a final state, use the original
    if final_state is None:
        logging.warning("All equilibration strategies failed, using original state")
        final_state = state_min
        S1 = final_state.entropy()
        S2 = S1
        improvement = 0

    logging.info("Improvement from equilibration: " + str(improvement))
    logging.info("Final entropy after equilibration: " + str(final_state.entropy()))

    logging.info("Saving blockstate...")
    block_state = final_state.get_bs()
    with open(snakemake.output.blocks, 'wb') as fh:
        dill.dump(block_state, fh, recurse=True)
    logging.info("Done!")