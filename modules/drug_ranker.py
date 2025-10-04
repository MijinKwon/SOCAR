import pandas as pd


def rank_drugs(cfg, pss_df):
    approved_df = pd.read_csv(cfg["data"]["approved_drugs"], sep="\t")
    approved_set = set(approved_df["DrugName"])


    pss_df = pss_df[pss_df["DrugName"].isin(approved_set)]
    ranked = pss_df.sort_values(by="PSS", ascending=False)
    ranked.to_csv(cfg["data"]["ranked_drugs"], sep="\t", index=False)