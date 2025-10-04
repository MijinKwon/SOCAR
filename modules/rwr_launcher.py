import pandas as pd
import numpy as np
import os
from utils.logger import get_logger

logger = get_logger("rwr_launcher")


def prepare_inputs(cfg, degs, network, drug_targets):
    """Prepare input files for RWR."""
    degs[["GeneSymbol", "Expression"]].to_csv(cfg["data"]["rwr_input"], sep="\t", index=False, header=True)
    network.to_csv(cfg["data"]["rwr_network"], sep="\t", index=False, header=True)
    drug_targets.to_csv(cfg["data"]["rwr_drug_targets"], sep="\t", index=False, header=True)


def build_transition(network_path: str):
    """Build row-stochastic transition matrix A and node index maps from a TSV network file."""
    network_df = pd.read_csv(network_path, sep="\t", header=0, names=["source", "target", "type", "interaction"])
    nodes = sorted(set(network_df["source"]) | set(network_df["target"]))
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}
    idx_to_node = {idx: node for node, idx in node_to_idx.items()}
    n = len(nodes)
    adj = np.zeros((n, n))
    for _, row in network_df.iterrows():
        i = node_to_idx[row["source"]]
        j = node_to_idx[row["target"]]
        adj[i, j] = 1
        # adj[j, i] = 1
    row_sums = adj.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1e-8
    A = adj / row_sums
    return A, node_to_idx, idx_to_node


def _rwr_from_seeds(A, seed_indices, alpha, epsilon, max_iter, seed_total_mass=30.0):
    n = A.shape[0]
    p0 = np.zeros(n)
    if not seed_indices:
        return p0
    p0[seed_indices] = 1.0
    # Normalize and scale so that total seed mass equals seed_total_mass
    p0 = (p0 / p0.sum()) * float(seed_total_mass)
    p = p0.copy()
    for _ in range(max_iter):
        p_new = (1 - alpha) * A @ p + alpha * p0
        if np.linalg.norm(p_new - p, 1) < epsilon:
            return p_new
        p = p_new
    return p


def run_per_drug_rwr(cfg, network_df: pd.DataFrame, drug_targets_df: pd.DataFrame):
    """Run RWR per drug using its targets as seeds and save results under data.rwr_per_drug_dir.

    Output: one TSV per drug: <dir>/<DrugName>.tsv with columns [GeneSymbol, Score]
    """
    # Build transition from in-memory DF to avoid re-reading
    nodes = sorted(set(network_df.iloc[:, 0]) | set(network_df.iloc[:, 1]))
    node_to_idx = {n: i for i, n in enumerate(nodes)}
    idx_to_node = {i: n for n, i in node_to_idx.items()}
    n = len(nodes)
    adj = np.zeros((n, n))
    for _, row in network_df.iterrows():
        s = row.iloc[0]
        t = row.iloc[1]
        if s in node_to_idx and t in node_to_idx:
            adj[node_to_idx[s], node_to_idx[t]] = 1
    row_sums = adj.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1e-8
    A = adj / row_sums

    alpha = cfg["rwr"]["restart_prob"]
    epsilon = cfg["rwr"]["convergence_threshold"]
    max_iter = cfg["rwr"].get("max_iter", 100)
    seed_total_mass = cfg["rwr"].get("seed_total_mass", 30)

    out_dir = cfg["data"].get("rwr_per_drug_dir", "output/random_walk_result")
    os.makedirs(out_dir, exist_ok=True)

    # Load approved drug list directly from config
    approved_path = cfg["data"]["approved_drugs"]
    approved_df = pd.read_csv(approved_path, sep="\t")
    approved_set = set(approved_df["DrugName"].astype(str).str.lower())

    # Map drug -> targets present in the network
    results = {}
    seeds_log = []  # to record actually mapped seeds per drug
    for _, r in drug_targets_df.iterrows():
        drug = str(r["DrugName"]).strip()
        targets = [t for t in str(r["TargetGene"]).split("|") if t]
        if drug.lower() not in approved_set:
            continue
        mapped_targets = [t for t in targets if t in node_to_idx]
        seed_idx = [node_to_idx[t] for t in mapped_targets]
        if not seed_idx:
            logger.warning(f"Skipping {drug}: no targets found on network")
            continue
        p = _rwr_from_seeds(A, seed_idx, alpha, epsilon, max_iter, seed_total_mass)
        df = pd.DataFrame({
            "GeneSymbol": [idx_to_node[i] for i in range(n)],
            "Score": p
        }).sort_values("Score", ascending=False)
        results[drug] = df
        df.to_csv(os.path.join(out_dir, f"{drug}.tsv"), sep="\t", index=False)
        seeds_log.append({
            "DrugName": drug,
            "SeedCount": len(mapped_targets),
            "Seeds": "|".join(mapped_targets),
        })
    logger.info(f"Per-drug RWR results saved to {out_dir}")
    # Write seeds used log
    seeds_path = cfg["data"].get("rwr_seeds_used", "processed_data/rwr_seeds_used.tsv")
    os.makedirs(os.path.dirname(seeds_path), exist_ok=True)
    pd.DataFrame(seeds_log).to_csv(seeds_path, sep="\t", index=False)
    logger.info(f"Seeds used per drug saved to {seeds_path}")
    return results
