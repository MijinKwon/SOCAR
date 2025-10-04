import os
from collections import deque
import pandas as pd

from utils.logger import get_logger


logger = get_logger("rmg_extractor")


def _normalize_network_cols(df: pd.DataFrame) -> pd.DataFrame:
    cols = {c.lower(): c for c in df.columns}
    source_col = cols.get("sourcegene") or cols.get("source") or list(df.columns)[0]
    target_col = cols.get("targetgene") or cols.get("target") or list(df.columns)[1]
    type_col = cols.get("pathwaytype") or cols.get("type") or (list(df.columns)[2] if len(df.columns) > 2 else None)
    out = pd.DataFrame({
        "source": df[source_col].astype(str),
        "target": df[target_col].astype(str),
    })
    if type_col is not None:
        out["type"] = df[type_col].astype(str)
    else:
        out["type"] = ""
    return out


def _largest_connected_component_undirected(nodes, edges):
    # Build undirected adjacency
    adj = {n: set() for n in nodes}
    for u, v in edges:
        if u in adj and v in adj:
            adj[u].add(v)
            adj[v].add(u)

    visited = set()
    best_comp = set()
    for start in nodes:
        if start in visited:
            continue
        # BFS
        comp = set()
        q = deque([start])
        visited.add(start)
        while q:
            x = q.popleft()
            comp.add(x)
            for y in adj[x]:
                if y not in visited:
                    visited.add(y)
                    q.append(y)
        if len(comp) > len(best_comp):
            best_comp = comp
    return best_comp


def extract_rmgs(cfg, degs_df: pd.DataFrame, network_df: pd.DataFrame, acgs_df: pd.DataFrame) -> set:
    """Extract Resistance Module Genes (RMGs) as LCC of induced subgraph on PRGs = DEGs ∪ ACGs.

    Saves one gene per line to cfg['data']['rmgs'] and returns the set of RMGs.
    """
    path = cfg["data"].get("rmgs", "processed_data/RMGs.txt")

    net = _normalize_network_cols(network_df)
    net = net[net["type"].str.lower().isin(["signaling", "regulatory"]) | (net["type"] == "")]

    deg_col = "GeneSymbol" if "GeneSymbol" in degs_df.columns else degs_df.columns[0]
    prg = set(degs_df[deg_col].astype(str).tolist())
    if acgs_df is not None and not acgs_df.empty:
        acg_col = "GeneSymbol" if "GeneSymbol" in acgs_df.columns else acgs_df.columns[0]
        prg |= set(acgs_df[acg_col].astype(str).tolist())

    # Induced nodes and edges
    nodes = set()
    edges = []
    for _, r in net.iterrows():
        s = r["source"]
        t = r["target"]
        if s in prg and t in prg:
            nodes.add(s)
            nodes.add(t)
            edges.append((s, t))

    if not nodes:
        logger.warning("No PRG nodes found in network; RMGs set is empty.")
        rmgs = set()
    else:
        lcc = _largest_connected_component_undirected(nodes, edges)
        rmgs = set(lcc)

    logger.info(f"Number of RMGs: {len(rmgs)}")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for g in sorted(rmgs):
            f.write(f"{g}\n")
    logger.info(f"RMGs saved to {path}")
    return rmgs

