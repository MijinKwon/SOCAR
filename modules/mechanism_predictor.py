import os
from typing import Set
import pandas as pd

from utils.logger import get_logger
from modules.enrichment import perform_go_enrichment


logger = get_logger("mechanism_predictor")


def predict_sensitization_mechanisms(
    cfg,
    ranked_df: pd.DataFrame,
    rmgs: Set[str],
    top_drugs: int = 5,
    top_rmgs: int = 50,
):
    """For top-ranked drugs, use saved per-drug RWR files to select top RMGs and run GO enrichment.

    - Reads per-drug RWR scores from cfg['data']['rwr_per_drug_dir']/<Drug>.tsv
    - Writes top-RMGs to cfg['data']['mechanisms_dir']/<Drug>_top_RMGs.tsv
    - Runs GO enrichment over those top RMGs using modules.enrichment.perform_go_enrichment,
      directing outputs under data.gsea_dir_mechanisms/<Drug>/
    """
    mech_dir = cfg["data"].get("mechanisms_dir", "output/mechanisms")
    os.makedirs(mech_dir, exist_ok=True)
    per_dir = cfg["data"].get("rwr_per_drug_dir", "output/random_walk_result")

    # Select top drugs from ranked_df
    drug_col = "DrugName" if "DrugName" in ranked_df.columns else ranked_df.columns[0]
    top = ranked_df.head(top_drugs)[drug_col].astype(str).tolist()

    for drug in top:
        drug_file = os.path.join(per_dir, f"{drug}.tsv")
        if not os.path.exists(drug_file):
            logger.warning(f"Per-drug RWR file not found for {drug}: {drug_file}")
            continue
        df = pd.read_csv(drug_file, sep="\t")
        df = df[df["GeneSymbol"].isin(rmgs)].sort_values("Score", ascending=False).head(top_rmgs)
        out_df = df[["GeneSymbol", "Score"]].copy()
        out_path = os.path.join(mech_dir, f"{drug}_top_RMGs.tsv")
        out_df.to_csv(out_path, sep="\t", index=False)
        logger.info(f"Mechanism RMGs for {drug} saved to {out_path}")

        # Run GO enrichment for top RMGs, placing outputs in a separate drug_mechanisms directory
        data_cfg = cfg.get("data", {})
        base_mech_dir = data_cfg.get(
            "gsea_dir_mechanisms",
            os.path.join(data_cfg.get("gsea_dir", "gene_set_enrichment_analysis"), "drug_mechanisms"),
        )
        per_drug_dir = os.path.join(base_mech_dir, drug)
        new_cfg = dict(cfg)
        new_cfg["data"] = dict(cfg.get("data", {}))
        # Set generic gsea_dir so enrichment writes directly into per-drug folder
        new_cfg["data"]["gsea_dir"] = per_drug_dir
        try:
            # Use a non-DEG/RMG description to avoid enrichment routing into DEGs/ or RMGs/
            perform_go_enrichment(out_df["GeneSymbol"].astype(str).tolist(), "MECHANISMS", new_cfg)
        except Exception as e:
            logger.warning(f"GO enrichment skipped for {drug}: {e}")
