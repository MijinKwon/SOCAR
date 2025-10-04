import pandas as pd
from sklearn.metrics import roc_auc_score


def evaluate(cfg):
    ranked = pd.read_csv(cfg["data"]["ranked_drugs"], sep="\t")
    gold = pd.read_csv(cfg["data"]["gold_standard"], sep="\t")
    approved = pd.read_csv(cfg["data"]["approved_drugs"], sep="\t")

    approved_set = set(approved["DrugName"])
    ranked = ranked[ranked["DrugName"].isin(approved_set)].copy()

    gold_set = set(gold["DrugName"])
    ranked["Label"] = ranked["DrugName"].apply(lambda x: 1 if x in gold_set else 0)

    auroc = roc_auc_score(ranked["Label"], ranked["PSS"])
    with open(cfg["data"]["evaluation_result"], "w") as f:
        f.write(f"AUROC: {auroc:.4f}\n")
