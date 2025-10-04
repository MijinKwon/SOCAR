from config import load_config
from utils.logger import get_logger
from utils.file_io import ensure_directories_exist
from modules.preprocessing import load_inputs
from modules.rwr_launcher import prepare_inputs, run_per_drug_rwr
from modules.pss_calculator import calculate_pss_from_per_drug
from modules.drug_ranker import rank_drugs
from modules.evaluator import evaluate
from modules.acg_extractor import extract_acgs
from modules.rmg_extractor import extract_rmgs
from modules.mechanism_predictor import predict_sensitization_mechanisms
from modules.enrichment import perform_go_enrichment


def main():
    cfg = load_config()
    logger = get_logger("main")
    ensure_directories_exist(["processed_data", "output", "logs"])

    # logger.info("Loading input files...")
    degs, network, drug_targets = load_inputs(cfg)

    # 1) Extract ACGs and RMGs (PRGs = DEGs ∪ ACGs; LCC := RMGs)
    # logger.info("Extracting ACGs and RMGs from DEGs and network...")
    acgs_df = extract_acgs(cfg, degs, network)
    rmgs_set = extract_rmgs(cfg, degs, network, acgs_df)

    # Optional: GO enrichment for DEGs and RMGs before RWR
    enr_cfg = cfg.get("analysis", {}).get("enrichment", {})
    if enr_cfg.get("enabled", False):
        logger.info("Running GO enrichment for DEGs and RMGs...")
        deg_col = "GeneSymbol" if "GeneSymbol" in degs.columns else degs.columns[0]
        perform_go_enrichment(degs[deg_col].astype(str).tolist(), "DEGs", cfg)
        perform_go_enrichment(sorted(rmgs_set), "RMGs", cfg)

    # 2) Prepare I/O artifacts and run per‑drug RWR (targets as seeds)
    logger.info("Preparing inputs and running per-drug RWR...")
    prepare_inputs(cfg, degs, network, drug_targets)
    run_per_drug_rwr(cfg, network, drug_targets)

    # 3) Compute PSS using selected gene set per config (DEGs or RMGs)
    logger.info("Calculating Personalized Sensitivity Scores (PSS)...")
    gene_option = (cfg.get("analysis", {}).get("gene_option", "RMGs")).upper()
    if gene_option in ("DEG", "DEGS"):
        deg_col = "GeneSymbol" if "GeneSymbol" in degs.columns else degs.columns[0]
        allowed_genes = set(degs[deg_col].astype(str).tolist())
    else:
        allowed_genes = set(rmgs_set)
    pss_df = calculate_pss_from_per_drug(cfg, allowed_genes, drug_targets)

    logger.info("Ranking drugs based on PSS scores...")
    rank_drugs(cfg, pss_df)

    # Optional: mechanism prediction via per-drug RWR
    analysis_cfg = cfg.get("analysis", {})
    if analysis_cfg.get("enable_mechanism_prediction", False):
        import pandas as pd
        logger.info("Predicting sensitization mechanisms for top-ranked drugs...")
        ranked_df = pd.read_csv(cfg["data"]["ranked_drugs"], sep="\t")
        top_drugs = int(analysis_cfg.get("top_drugs", 5))
        top_rmgs = int(analysis_cfg.get("top_rmgs", 50))
        predict_sensitization_mechanisms(
            cfg,
            ranked_df,
            rmgs_set,
            top_drugs=top_drugs,
            top_rmgs=top_rmgs,
        )

    logger.info("Evaluating results using AUROC against gold standard...")
    evaluate(cfg)

    logger.info("Done.")


if __name__ == "__main__":
    main()
