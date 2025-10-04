# Input File Schema

This document describes the format of each input file in the `input_data/` directory.

---

## DEGs.txt

| Column       | Description                        |
|--------------|------------------------------------|
| GeneSymbol   | HGNC-approved gene symbol          |
| Expression   | DEG_Up or DEG_Down                 |

- List of differentially expressed genes.
- Used to extract Activity Change Genes (ACGs).

---

## drug_target_associations.txt

| Column     | Description                                          |
|------------|------------------------------------------------------|
| DrugName   | Name of the drug                                     |
| TargetGene | One or more gene symbols, separated by pipe    |

- Known drug–target associations.
- Multiple targets for one drug are allowed.

---

## approved_drugs.txt

| Column     | Description                |
|------------|----------------------------|
| DrugName   | FDA-approved drug name     |

- FDA-approved drug list used to filter candidates.

---

## gold_standard_drugs.txt

| Column     | Description                          |
|------------|--------------------------------------|
| DrugName   | Known sensitizer drug name           |
| TargetGene | Known gene target (HGNC symbol)      |

- Used for AUROC performance benchmarking.

---

## regulatory_signaling_network.txt

| Column        | Description                                                      |
|---------------|------------------------------------------------------------------|
| SourceGene    | Upstream entity (HGNC gene symbol)                               |
| TargetGene    | Downstream gene (HGNC symbol)                                    |
| PathwayType   | `signaling` or `regulatory`                                      |
| Interaction   | One of: activate, inhibit, express, repress, etc.                |

- Includes curated interactions with FDR < 0.05.
- Constructed using prior KGML extraction pipelines and publicly available databases.
