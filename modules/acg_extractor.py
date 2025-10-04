import os
import pandas as pd

from utils.logger import get_logger


logger = get_logger("acg_extractor")


def _normalize_network_cols(df: pd.DataFrame) -> pd.DataFrame:
    cols = {c.lower(): c for c in df.columns}
    # Prefer common header names; fallback to positional if needed
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


def extract_acgs(cfg, degs_df: pd.DataFrame, network_df: pd.DataFrame) -> pd.DataFrame:
    """Extract Activity Change Genes (ACGs) as non-DEG signaling neighbors of DEGs.

    Writes a two-column TSV to cfg['data']['acgs'] with [GeneSymbol, Count].
    Returns the ACG DataFrame.
    """
    path = cfg["data"].get("acgs", "processed_data/ACGs.txt")

    net = _normalize_network_cols(network_df)
    net_sig = net[net["type"].str.lower().eq("signaling")]
    if net_sig.empty:
        # If type column is missing or not labeled, fallback to all edges
        logger.warning("No 'signaling' edges found; using all network edges for ACG extraction.")
        net_sig = net

    deg_col = "GeneSymbol" if "GeneSymbol" in degs_df.columns else degs_df.columns[0]
    degs = set(degs_df[deg_col].astype(str).tolist())

    # Build adjacency: out-neighbors per source
    adj = {}
    for _, r in net_sig.iterrows():
        s = r["source"]
        t = r["target"]
        adj.setdefault(s, set()).add(t)

    acg_counts = {}
    matched = 0
    for g in degs:
        nbrs = adj.get(g)
        if not nbrs:
            continue
        matched += 1
        for n in nbrs:
            if n in degs:
                continue
            acg_counts[n] = acg_counts.get(n, 0) + 1

    logger.info(f"DEGs with signaling neighbors: {matched}")
    logger.info(f"Number of ACGs identified: {len(acg_counts)}")

    acg_df = pd.DataFrame({"GeneSymbol": list(acg_counts.keys()), "Count": list(acg_counts.values())})
    acg_df = acg_df.sort_values(["GeneSymbol"]).reset_index(drop=True)

    os.makedirs(os.path.dirname(path), exist_ok=True)
    acg_df.to_csv(path, sep="\t", index=False)
    logger.info(f"ACGs saved to {path}")
    return acg_df

