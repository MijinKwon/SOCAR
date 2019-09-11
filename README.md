# SOCAR
Discovery of sensitizers to overcome cancer drug resistance

[Environment]
python-3.7

[How to use]
1. INPUT
   - The following files are required as input files in the 'input_data' directory.
     - approved_drugs.txt : the list of FDA-approved drugs
     - drug_target_association.txt : physical interactions between drugs and target proteins
     - regulatory_signaling_network : physical and directed molecular interactions
     - DEGs.txt : differentially expressed genes between cancer cells in sensitive or resistant to a drug
     - gold_standard_drugs: (for performance evaluation) known sensitizers with sensitization effects for a resistant drug

2. MODEL
   - Run 'SOCAR.py' code.

3. OUTPUT
   - PSSs of test drugs are calculated in the 'model_evaluation' directory.
   - Model performance is evaluated using AUROC measurment in the 'model_evaluation' directory.
   - Sensitization mechanisms of candidate sensitizers are predicted in the 'gene_set_enrichment_analysis' directory.
