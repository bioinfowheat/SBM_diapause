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
from multiprocessing import cpu_count

from graph_tool.all import *
import numpy as np

def run_parallel_tempering_mcmc(state_min, niter, n_chains=8):
    """
    Run parallel tempering MCMC for better multi-core utilization
    """
    logging.info(f"Running parallel tempering MCMC with {n_chains} chains...")
    
    # Create temperature schedule
    betas = np.linspace(1.0, 0.1, n_chains)
    logging.info(f"Temperature schedule (betas): {betas}")
    
    # Create multiple states at different temperatures
    states = [state_min.copy() for _ in range(n_chains)]
    
    # Storage for collecting partitions
    bs = []
    h = [np.zeros(state_min.g.num_vertices() + 1) for s in state_min.get_levels()]
    
    def collect_partitions_pt(states):
        # Collect from the ground state (beta=1.0, first state)
        s = states[0]
        bs.append(s.get_bs())
        for l, sl in enumerate(s.get_levels()):
            B = sl.get_nonempty_B()
            h[l][B] += 1
    
    try:
        # Try parallel tempering
        S1 = state_min.entropy()
        multicanonical_equilibrate(states, betas, 
                                  force_niter=niter,
                                  mcmc_args=dict(niter=10),
                                  callback=collect_partitions_pt, 
                                  verbose=True)
        
        # Get the ground state result
        final_state = states[0]
        S2 = final_state.entropy()
        
        return bs, h, S1, S2, final_state
        
    except Exception as e:
        logging.warning(f"Parallel tempering failed: {e}")
        logging.info("Falling back to sequential MCMC...")
        return None

def run_sequential_mcmc(state_min, niter):
    """
    Fallback to sequential MCMC
    """
    logging.info("Running sequential MCMC...")
    
    bs = []
    h = [np.zeros(state_min.g.num_vertices() + 1) for s in state_min.get_levels()]

    def collect_partitions(s):
        bs.append(s.get_bs())
        for l, sl in enumerate(s.get_levels()):
            B = sl.get_nonempty_B()
            h[l][B] += 1

    S1 = state_min.entropy()
    mcmc_equilibrate(state_min, force_niter=niter, mcmc_args=dict(niter=10),
                    callback=collect_partitions, verbose=True)
    
    S2 = state_min.entropy()
    
    return bs, h, S1, S2, state_min

def run_multiple_independent_chains(state_min, niter, n_chains=4):
    """
    Alternative: Run multiple independent MCMC chains and combine results
    """
    from multiprocessing import Pool
    import copy
    
    logging.info(f"Running {n_chains} independent MCMC chains in parallel...")
    
    def run_single_chain(args):
        state, chain_id, chain_niter = args
        
        # Set up logging for this chain
        chain_logger = logging.getLogger(f'chain_{chain_id}')
        
        bs_chain = []
        h_chain = [np.zeros(state.g.num_vertices() + 1) for s in state.get_levels()]
        
        def collect_partitions_chain(s):
            bs_chain.append(s.get_bs())
            for l, sl in enumerate(s.get_levels()):
                B = sl.get_nonempty_B()
                h_chain[l][B] += 1
        
        try:
            mcmc_equilibrate(state, force_niter=chain_niter, mcmc_args=dict(niter=10),
                           callback=collect_partitions_chain, verbose=False)
            
            return bs_chain, h_chain, state.entropy()
        except Exception as e:
            chain_logger.error(f"Chain {chain_id} failed: {e}")
            return [], h_chain, float('inf')
    
    # Prepare arguments for parallel execution
    chain_niter = max(1, niter // n_chains)
    args = [(state_min.copy(), i, chain_niter) for i in range(n_chains)]
    
    try:
        # Run chains in parallel
        with Pool(n_chains) as pool:
            results = pool.map(run_single_chain, args)
        
        # Combine results
        bs_combined = []
        h_combined = [np.zeros(state_min.g.num_vertices() + 1) for s in state_min.get_levels()]
        entropies = []
        
        for bs_chain, h_chain, entropy in results:
            bs_combined.extend(bs_chain)
            for l in range(len(h_combined)):
                h_combined[l] += h_chain[l]
            entropies.append(entropy)
        
        # Find best final state
        best_entropy = min(entropies)
        logging.info(f"Best entropy from parallel chains: {best_entropy}")
        
        return bs_combined, h_combined, state_min.entropy(), best_entropy, state_min
        
    except Exception as e:
        logging.error(f"Parallel chain execution failed: {e}")
        return None

if __name__ == '__main__':
    
    # Check system capabilities
    logging.info(f"Graph-tool OpenMP enabled: {openmp_enabled()}")
    logging.info(f"Available CPU cores: {cpu_count()}")
    
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

    # Determine strategy based on graph size and available cores
    n_vertices = g.num_vertices()
    n_cores = min(8, cpu_count())
    
    logging.info("Starting MCMC...")
    
    # Strategy 1: Try parallel tempering for larger graphs
    if n_vertices > 1000 and openmp_enabled():
        result = run_parallel_tempering_mcmc(state_min, snakemake.params.niter, n_chains=n_cores)
        
        if result is not None:
            bs, h, S1, S2, final_state = result
            logging.info("Successfully used parallel tempering MCMC")
        else:
            # Fallback to multiple independent chains
            logging.info("Trying multiple independent chains...")
            result = run_multiple_independent_chains(state_min, snakemake.params.niter, n_chains=4)
            
            if result is not None:
                bs, h, S1, S2, final_state = result
                logging.info("Successfully used multiple independent chains")
            else:
                # Final fallback to sequential
                bs, h, S1, S2, final_state = run_sequential_mcmc(state_min, snakemake.params.niter)
                logging.info("Used sequential MCMC as fallback")
    
    # Strategy 2: For medium graphs, try multiple independent chains
    elif n_vertices > 500:
        result = run_multiple_independent_chains(state_min, snakemake.params.niter, n_chains=4)
        
        if result is not None:
            bs, h, S1, S2, final_state = result
            logging.info("Successfully used multiple independent chains")
        else:
            bs, h, S1, S2, final_state = run_sequential_mcmc(state_min, snakemake.params.niter)
            logging.info("Used sequential MCMC as fallback")
    
    # Strategy 3: For small graphs, use sequential (most efficient)
    else:
        bs, h, S1, S2, final_state = run_sequential_mcmc(state_min, snakemake.params.niter)
        logging.info("Used sequential MCMC (optimal for small graphs)")

    # Post-processing (same as original)
    logging.info("Processing results...")
    pmode = PartitionModeState(bs, nested=True, converge=True)
    pv = pmode.get_marginal(g)

    # Get consensus estimate
    bs_max = pmode.get_max_nested()
    state = final_state.copy(bs=bs_max)
    S2_final = state.entropy()

    logging.info("Description length improvement in MCMC: " + str(S2_final - S1))
    logging.info("Final entropy after MCMC: " + str(state.entropy()))

    logging.info("Saving blockstate...")
    block_state = state.get_bs()
    with open(snakemake.output.blocks, 'wb') as fh:
        dill.dump(block_state, fh, recurse=True)

    logging.info("Saving block number histogram...")
    with open(snakemake.output.hist, 'wb') as fh:
        dill.dump(h, fh, recurse=True)

    logging.info("Done!")