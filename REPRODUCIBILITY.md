# Reproducibility Guide for SOCAR

This document provides a step-by-step guide to reproduce the main results of the manuscript:

“SOCAR: Network-based computational framework to overcome acquired tamoxifen resistance of MCF7 cells.”

---

## 1. Environment Setup

### Requirements
- Python ≥ 3.10 (tested on 3.13.2)
- Packages: pandas, numpy, PyYAML, scikit-learn
- Optional (for enrichment): gseapy

### Install
```
pip install -r requirements.txt
```
Or manually:
```
pip install pandas numpy pyyaml scikit-learn gseapy
```

---

## 2. Input Data Preparation

All input formats and sources are described in:
- [input_data/INPUT_SCHEMA.md](input_data/INPUT_SCHEMA.md)
- [input_data/DATA_SOURCES.md](input_data/DATA_SOURCES.md)


### Minimal Example Files
Sample inputs (for quick testing) can be organized as:
```
input_data/
├── DEGs_sample.txt
├── drug_target_associations_sample.txt
├── regulatory_signaling_network_sample.txt
├── approved_drugs.txt
└── gold_standard_drugs.txt
```
To reproduce the full-scale analysis (MCF7 dataset), request access to the authors.

---

## 3. Configuration

Key parameters are defined in `config.yaml`.

Default important settings:
```
rwr:
  restart_prob: 0.7
  convergence_threshold: 1e-5
  max_iter: 100
  seed_total_mass: 30
analysis:
  gene_option: RMGs
```

---

## 4. Run the Pipeline

```
python main.py
```

Expected runtime: ~2–5 min on standard CPU for sample data.

---

## 5. Outputs

Main results are written under `output/`:
```
output/
├── pss_scores.txt          # Drug-wise potential sensitivity scores
├── drug_rankings.txt       # Ranked drugs by PSS
├── auroc_result.txt        # Evaluation result (expected AUROC ≈ 0.93)
├── random_walk_result/     # Per-drug RWR profiles
└── mechanisms/             # Mechanism-level top-RMG results
```
Intermediate data under `processed_data/` are automatically generated.

---

## 6. Verification

You can check the performance with:
```
cat output/auroc_result.txt
```
Expected result (based on MCF7/TAMR dataset):
```
AUROC = 0.93 ± 0.01
```

---

## 7. Optional: GO Enrichment Analysis

If `gseapy` is installed, the following directories will be generated:
```
gene_set_enrichment_analysis/
├── DEGs/
├── RMGs/
└── drug_mechanisms/
```
Each contains GO terms enriched for the respective gene sets.

---

## 8. Citation

If you use this framework, please cite:

> Kwon, M. (2025). SOCAR: Network-based computational framework to overcome acquired tamoxifen resistance of MCF7 cells.  
> GitHub: [https://github.com/MijinKwon/SOCAR](https://github.com/MijinKwon/SOCAR)  
> License: CC-BY-4.0


---

## 9. Contact

For questions or access to the full input datasets, please contact:

Mijin Kwon — via GitHub Issues

