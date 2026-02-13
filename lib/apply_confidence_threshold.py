#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import sys

# ============================================================
# USER SETTINGS
# ============================================================

DEMUX_DIR = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/ZC-Basecall/demux"   # <-- EDIT
INPUT_TSV = os.path.join(DEMUX_DIR, "demux_per-read_modifications.tsv")
OUTPUT_CSV = os.path.join(DEMUX_DIR, "overall_demux_results_confident.csv")

CONF_THRESHOLD = 0.98   # e.g. 0.8, 0.9, 0.95

# ============================================================
# Helper functions
# ============================================================

def parse_class1_prob(prob_string):
    """
    Extract P(class1) from 'class_probs'
    Assumes format: 'P0,P1'
    """
    try:
        return float(prob_string.split(",")[1])
    except Exception:
        return np.nan


def apply_confidence_filter(df, threshold):
    """
    Keep reads that meet confidence requirements
    """
    df = df.copy()
    df["p_class1"] = df["class_probs"].apply(parse_class1_prob)

    confident = (
        ((df["class_pred"] == 1) & (df["p_class1"] >= threshold)) |
        ((df["class_pred"] == 0) & (df["p_class1"] <= (1 - threshold)))
    )

    return df[confident]


def calculate_overall_results(df):
    """
    Aggregate confident reads
    """
    results = (
        df.groupby(["sample_id", "barcode_pair"])
          .agg(
              Total_Reads=("read_id", "size"),
              Number_of_1s=("class_pred", lambda x: (x == 1).sum()),
              Number_of_0s=("class_pred", lambda x: (x == 0).sum())
          )
          .reset_index()
    )

    results["Percentage_1"] = 100 * results["Number_of_1s"] / results["Total_Reads"]
    results["Percentage_0"] = 100 * results["Number_of_0s"] / results["Total_Reads"]

    return results


# ============================================================
# Main
# ============================================================

print(f"[INFO] Loading demuxed per-read file: {INPUT_TSV}")
df = pd.read_csv(INPUT_TSV, sep="\t")

print(f"[INFO] Applying confidence threshold: {CONF_THRESHOLD}")
df_conf = apply_confidence_filter(df, CONF_THRESHOLD)

print(f"[INFO] Reads before filter: {len(df)}")
print(f"[INFO] Reads after filter:  {len(df_conf)}")

overall_results = calculate_overall_results(df_conf)

overall_results.to_csv(OUTPUT_CSV, index=False)

print(f"[DONE] Confident overall results written to:")
print(f"       {OUTPUT_CSV}")

