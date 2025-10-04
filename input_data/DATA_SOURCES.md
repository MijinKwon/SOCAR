# Data Sources

This document summarizes the origin and curation strategy of all input files in the SOCAR framework.

---

## DEGs.txt

- **Source**: GEO dataset [GSE67916](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67916)
- **Samples**: 8 tamoxifen-sensitive MCF7 lines, 3 tamoxifen-resistant MCF7/TAMR lines
- **Platform**: Affymetrix Human Genome U133 Plus 2.0 (GPL570)
- **Processing**: GCRMA normalization and LIMMA for differential expression analysis
- **Note**: Only the DEG result file is included (`DEGs.txt`). Raw CEL files are **not** included in this repository.

---

## regulatory_signaling_network.txt

- **Original Network**: Constructed in prior study [60] by the authors
- **Integrated Databases**:
  - **KEGG**: Extracted via KEGGgraph from KGML (activation, inhibition, etc.)
  - **SIGNOR**: Signaling and transcriptional regulation mechanisms
  - **TRANSFAC**: Transcriptional regulatory edges
- **Curation**:
  - Includes only directed human protein–gene interactions
  - All gene symbols converted to HGNC-approved names

---

## drug_target_associations.txt

- **Source**: DrugBank database
- **Content**: Known drug–target gene mappings
- **Format**: Multiple targets per drug allowed (pipe-separated)

---

## approved_drugs.txt

- **Source**: DrugBank FDA-approved drugs list
- **Content**: Filter of drugs approved by FDA at time of analysis

---

## gold_standard_drugs.txt

- **Source**: Literature-curated sensitizer drugs associated with tamoxifen resistance
- **Usage**: Used for validation (AUROC) of SOCAR predictions
- **Note**: 13 sensitizers were compiled; 6 have DrugBank therapeutic indications, and 4 are FDA-approved
