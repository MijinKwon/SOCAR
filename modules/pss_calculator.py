import os
import glob
import pandas as pd


def calculate_pss_from_per_drug(cfg, allowed_genes_set, drug_targets_df):
    """Calculate PSS using per-drug RWR results restricted to a given gene set.

    For each drug file, sum scores over genes in `allowed_genes_set` and normalize by
    the number of actually mapped target seeds for that drug (from rwr_seeds_used); falls back
    to declared target count if seed log is missing. Writes [DrugName, PSS].
    """
    per_dir = cfg["data"].get("rwr_per_drug_dir", "output/random_walk_result")
    out_path = cfg["data"]["pss_output"]

    # Load actually mapped seed counts per drug (preferred denominator)
    seeds_path = cfg["data"].get("rwr_seeds_used", "processed_data/rwr_seeds_used.tsv")
    seed_counts = {}
    if os.path.exists(seeds_path):
        try:
            seeds_df = pd.read_csv(seeds_path, sep="\t")
            if "DrugName" in seeds_df.columns and "SeedCount" in seeds_df.columns:
                seed_counts = {str(r["DrugName"]): int(r["SeedCount"]) for _, r in seeds_df.iterrows()}
        except Exception:
            seed_counts = {}
    # Fallback: Build drug -> targets list from input (if seed_counts missing)
    tgt_map = {}
    if not seed_counts:
        for _, r in drug_targets_df.iterrows():
            drug = str(r["DrugName"]).strip()
            tgts = [t for t in str(r["TargetGene"]).split("|") if t]
            tgt_map[drug] = tgts

    # Determine valid drugs list based on available metadata
    valid_drugs = set(seed_counts.keys()) if seed_counts else set(tgt_map.keys()) if tgt_map else None

    rows = []
    for path in glob.glob(os.path.join(per_dir, "*.tsv")):
        drug = os.path.basename(path).replace(".tsv", "")
        if valid_drugs is not None and drug not in valid_drugs:
            continue
        df = pd.read_csv(path, sep="\t")
        # sum over allowed genes (e.g., DEGs or RMGs)
        df = df[df["GeneSymbol"].isin(allowed_genes_set)]
        s = float(df["Score"].sum()) if not df.empty else 0.0
        # denom = None
        # if drug in seed_counts:
        #     denom = float(seed_counts[drug]) if seed_counts[drug] > 0 else 1.0
        # else:
        #     # fallback to raw target count if no seed log is available
        #     if not tgt_map:
        #         denom = 1.0
        #     else:
        #         denom = float(len(tgt_map.get(drug, []))) or 1.0
        # pss = s / denom
        pss = s
        rows.append({"DrugName": drug, "PSS": pss})

    pss_df = pd.DataFrame(rows)
    pss_df.to_csv(out_path, sep="\t", index=False)
    return pss_df
