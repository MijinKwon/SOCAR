import pandas as pd

def convert_sep_to_tab(filepath):
    with open(filepath, 'r') as f:
        data = f.read()

    data = data.replace('\\t', '\t')
    with open(filepath, 'w') as f:
        f.write(data)

def load_inputs(cfg):
    convert_sep_to_tab(cfg["data"]["degs"])
    convert_sep_to_tab(cfg["data"]["network"])
    convert_sep_to_tab(cfg["data"]["drug_targets"])

    degs = pd.read_csv(cfg["data"]["degs"], sep="\t")
    network = pd.read_csv(cfg["data"]["network"], sep="\t")
    drug_targets = pd.read_csv(cfg["data"]["drug_targets"], sep="\t")
    return degs, network, drug_targets
