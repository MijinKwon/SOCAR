import os
from typing import Iterable

from utils.logger import get_logger


logger = get_logger("enrichment")


def perform_go_enrichment(genes: Iterable[str], description: str, cfg) -> None:
    """Run Enrichr GO enrichment (if gseapy installed). Writes results to cfg['data']['gsea_dir'].

    This mirrors the legacy behavior but is optional; if gseapy is not available,
    it logs a warning and skips without failing the pipeline.
    """
    genes = [str(g) for g in genes if g]
    if not genes:
        logger.warning(f"No genes provided for enrichment: {description}")
        return

    try:
        import gseapy as gp  # optional dependency
    except Exception:
        logger.warning("gseapy not installed; skipping enrichment.")
        return

    data_cfg = cfg.get("data", {})
    base_dir = data_cfg.get("gsea_dir", "gene_set_enrichment_analysis")
    desc_upper = str(description).upper()
    if desc_upper.startswith("DEG"):
        outdir = data_cfg.get("gsea_dir_degs", os.path.join(base_dir, "DEGs"))
    elif desc_upper.startswith("RMG"):
        outdir = data_cfg.get("gsea_dir_rmgs", os.path.join(base_dir, "RMGs"))
    else:
        outdir = base_dir
    os.makedirs(outdir, exist_ok=True)
    gene_set = (
        cfg.get("analysis", {})
        .get("enrichment", {})
        .get("gene_set", "GO_Biological_Process_2015")
    )

    try:
        gp.enrichr(
            gene_list=genes,
            gene_sets=gene_set,
            organism="Human",
            outdir=outdir,
            cutoff=0.05,
        )
        logger.info(f"Enrichment for {description} written to {outdir}")
    except Exception as e:
        logger.warning(f"Enrichment failed for {description}: {e}")
