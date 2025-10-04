# SOCAR (Sensitizers to Overcome CAncer drug Resistance)
SOCAR is a reproducible Python pipeline that identifies potential sensitizer drugs by integrating transcriptomic alterations with a curated molecular interaction network.  
It implements the framework described in the manuscript *“SOCAR: Network-based computational framework to overcome acquired tamoxifen resistance of MCF7 cells.”*

## Manuscript
SOCAR: Network-based computational framework to overcome acquired tamoxifen resistance of MCF7 cells  

## Repository Structure
```
SOCAR/
├── main.py                         # Pipeline entry point
├── config.yaml                     # I/O paths and RWR hyperparameters
├── modules/                        # Pipeline modules
│   ├── preprocessing.py            # Input loading and normalization
│   ├── rwr_launcher.py             # Prepare inputs and run Python RWR
│   ├── pss_calculator.py           # Compute Potential Sensitivity Scores
│   ├── drug_ranker.py              # Rank drugs (filters to approved drugs)
│   ├── evaluator.py                # AUROC-based evaluation
│   ├── acg_extractor.py            # Extract ACGs (signaling neighbors of DEGs)
│   ├── rmg_extractor.py            # Extract RMGs (LCC of DEG∪ACG subgraph)
│   └── mechanism_predictor.py      # Per-drug RWR to suggest top RMGs
├── utils/
│   ├── file_io.py                  # Directory creation utilities
│   └── logger.py                   # Simple logger setup
├── input_data/                     # Input data
│   ├── INPUT_SCHEMA.md                 # Input file formats
│   ├── DATA_SOURCES.md                 # Data sources
│   ├── DEGs.txt
│   ├── approved_drugs.txt
│   ├── drug_target_associations.txt
│   ├── gold_standard_drugs.txt
│   └── regulatory_signaling_network.txt
├── processed_data/                 # Intermediate outputs (auto-generated)
│   ├── ACGs.txt
│   ├── RMGs.txt
│   ├── rwr_input.txt
│   ├── rwr_network.txt
│   └── rwr_drug_targets.txt
├── output/                         # Final outputs (auto-generated)
│   ├── pss_scores.txt
│   ├── drug_rankings.txt
│   ├── auroc_result.txt
│   ├── random_walk_result/         # Per-drug RWR scores
│   │   └── <DrugName>.tsv
│   └── mechanisms/
│       └── <DrugName>_top_RMGs.tsv
└── README.md
```
Note: RWR is implemented in Python (`modules/rwr_launcher.py`)

## Requirements
- Python: 3.10–3.13 (tested on 3.13.2)
- Packages: pandas, numpy, PyYAML, scikit-learn
- Optional: gseapy (for GO enrichment)

### Installation
```
pip install pandas numpy pyyaml scikit-learn gseapy
```

## Configuration (`config.yaml`)
- `data`: Input/output file paths
  - Inputs: `degs`, `network`, `drug_targets`, `approved_drugs`, `gold_standard`
  - Intermediate/outputs: `acgs`, `rmgs`, `rwr_input`, `rwr_network`, `rwr_drug_targets`, `pss_output`, `ranked_drugs`, `evaluation_result`, `rwr_per_drug_dir`, `rwr_seeds_used`, `mechanisms_dir`
- `rwr`: Random Walk hyperparameters
  - `restart_prob`: Restart probability (alpha)
  - `convergence_threshold`: Convergence threshold (L1 delta)
  - `max_iter`: Maximum iterations
  - `seed_total_mass`: Total mass assigned to initial seeds (default 30)
- `analysis`: Optional behavior flags
  - `gene_option`: `RMGs` (restrict PSS to RMGs) (default) or `DEGs`
  - `enable_mechanism_prediction`: Run per‑drug RWR for top drugs
  - `top_drugs`, `top_rmgs`: Counts for mechanism summaries
  - Enrichment outputs:
    - `data.gsea_dir`: base enrichment folder
    - `data.gsea_dir_degs`: DEGs enrichment folder (default `<gsea_dir>/DEGs`)
    - `data.gsea_dir_rmgs`: RMGs enrichment folder (default `<gsea_dir>/RMGs`)
    - `data.gsea_dir_mechanisms`: per‑drug mechanisms enrichment folder (default `<gsea_dir>/drug_mechanisms`)
Adjust paths and parameters to fit your environment.

## How to Run
Place input files under `input_data/` and run:
```
python main.py
```
During execution, `processed_data/`, `output/`, and `logs/` are created automatically and stepwise logs are printed.

## Pipeline Stages
- Input loading: `modules.preprocessing.load_inputs`
  - Converts literal `\t` to real tabs and loads TXTs.
- ACG extraction (optional): `modules.acg_extractor.extract_acgs`
  - Counts signaling neighbors of DEGs that are not DEGs.
- RMG extraction (optional): `modules.rmg_extractor.extract_rmgs`
  - Builds DEG∪ACG induced subgraph; uses largest connected component as RMGs.
- RWR prep/run: `modules.rwr_launcher.prepare_inputs`, `modules.rwr_launcher.run_per_drug_rwr`
  - Builds a transition matrix and, for each approved drug, performs Random Walk with Restart using its mapped target genes as seeds.
- PSS computation: `modules.pss_calculator.calculate_pss_from_per_drug`
  - For each drug: sum per‑drug RWR scores over the selected geneset (DEGs or RMGs)
- Drug ranking: `modules.drug_ranker.rank_drugs`
  - Filters to approved drugs and sorts by PSS (descending).
- Evaluation: `modules.evaluator.evaluate`
  - Computes AUROC against the gold-standard set.
- Mechanism prediction (optional): `modules.mechanism_predictor.predict_sensitization_mechanisms`
  - For top drugs, runs per‑drug RWR and reports top‑scoring RMGs per drug.

## Outputs
- `processed_data/ACGs.txt`: ACG list with counts
- `processed_data/RMGs.txt`: RMG list (one per line)
- `processed_data/rwr_seeds_used.tsv`: Mapped target seeds per drug (records the mapped target seeds per drug; used as reference for RWR seed mapping)
- `output/random_walk_result/<DrugName>.tsv`: Per-drug RWR scores
- `output/pss_scores.txt`: Drug-wise PSS
- `output/drug_rankings.txt`: Final ranking for approved drugs
- `output/auroc_result.txt`: AUROC score
- `output/mechanisms/<DrugName>_top_RMGs.tsv`: Top RMGs per top-ranked drug
- `gene_set_enrichment_analysis/DEGs/`: GO enrichment for DEGs
- `gene_set_enrichment_analysis/RMGs/`: GO enrichment for RMGs (global)
- `gene_set_enrichment_analysis/drug_mechanisms/<DrugName>/`: GO enrichment for top-RMGs per drug

## Methodology Details

- Data loading: reads tab-delimited inputs; literal "\t" are converted to real tabs. Required columns are described in `input_schema.md`.
- ACG extraction: signaling neighbors of DEGs not in the DEG set are counted; outputs `processed_data/ACGs.txt` with columns `GeneSymbol`, `Count`.
- RMG extraction: constructs the subgraph induced by PRGs = DEGs ∪ ACGs on signaling+regulatory edges and takes the largest connected component (undirected) as the resistance module; outputs `processed_data/RMGs.txt`.
- Per-drug RWR: for each approved drug, uses its target genes that exist in the network as seeds and iteratively updates probabilities until convergence.
  - Transition: directed adjacency row-normalized to a stochastic matrix `A`, where each row represents the outgoing probability distribution from a gene to its downstream interactors in the signaling/regulatory network.
  - Seeds: initial vector `p0` assigns equal mass to each mapped target, scaled so the total mass equals `rwr.seed_total_mass` (default 30).
  - Restart: at each iteration, update by `p_new = (1 - alpha) * A @ p + alpha * p0` with `alpha = rwr.restart_prob` (default 0.7), until convergence (`L1 < rwr.convergence_threshold`).
  - Artifacts: per-drug score files under `output/random_walk_result/` and the per-drug seed counts in `processed_data/rwr_seeds_used.tsv`.
- PSS computation: for each drug, sum per-drug RWR scores over the selected geneset (from `rwr_seeds_used.tsv`). The allowed geneset is controlled by `analysis.gene_option`:
  - `DEGs`: uses the DEG list
  - `RMGs`: uses the resistance module gene set
- Ranking: filters to approved drugs and sorts by descending PSS; writes `output/drug_rankings.txt`.
- Evaluation: labels drugs as positive if present in `gold_standard_drugs.txt` and computes AUROC on approved drugs only; writes `output/auroc_result.txt`. The default configuration reproduces the manuscript result with AUROC ≈ 0.93 on the Tamoxifen-resistant MCF7 dataset.
- Mechanism analysis: for top-ranked drugs, selects top `analysis.top_rmgs` RMGs by per-drug RWR score, writes `output/mechanisms/<Drug>_top_RMGs.tsv`, and runs GO enrichment into `gene_set_enrichment_analysis/drug_mechanisms/<Drug>/`.

Notes
- If `gseapy` is not installed, enrichment steps are skipped with a log message.
- Nodes with no outgoing edges are handled by a tiny epsilon in row sums to avoid division-by-zero.
- All inputs must use HGNC-approved gene symbols and be tab-delimited.

## Input Formats and Data Sources
- Schema: `input_schema.md`
- Sources: `DATA_SOURCES.md`
- Note: All inputs are tab-delimited text files (literal `\t` or actual tabs). Use HGNC-approved gene symbols.

### Input Data Availability

Due to database redistribution policies, only sample input files are included in this repository.  
Full datasets are publicly available and can be downloaded as described in `input_data/DATA_SOURCES.md`.  
Place the full versions under `input_data/` with the same filenames to reproduce the full results.

## License
This code is provided for academic research use only.  
For commercial or redistribution inquiries, please contact the corresponding author.

