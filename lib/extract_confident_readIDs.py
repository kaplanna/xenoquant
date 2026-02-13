#!/usr/bin/env python3

"""
extract_high_confidence_read_ids.py

Extract read IDs from a per-read modifications file where:
- class_pred == 1
- p_class1 >= PROB_THRESHOLD
- barcode_pair is in a specified allow-list

Outputs a plain text file with one read_id per line.
"""

import pandas as pd
import os

# ======================================================
# ================= USER SETTINGS ======================
# ======================================================


PROB_THRESHOLD = 0.99

INPUT_FILE = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251215_8L-Ext_Basecall_NB24R/PG-75-75-6L/demux/demux_per-read_modifications.tsv"   # CSV or TSV
OUTPUT_TXT = f'/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251215_8L-Ext_Basecall_NB24R/PG_Highconf_reads_{PROB_THRESHOLD}.txt'

BARCODE_PAIRS = [
    "None_FWD_NB24_REV",
]

FILE_SEPARATOR = "\t"

# ======================================================
# ===================== SCRIPT =========================
# ======================================================

def extract_p_class1(class_probs):
    """
    class_probs format: "p0,p1"
    returns p1 as float
    """
    try:
        return float(class_probs.split(",")[1])
    except Exception:
        return float("nan")

def main():
    print("[INFO] Loading per-read modifications file...")
    df = pd.read_csv(INPUT_FILE, sep=FILE_SEPARATOR)

    required_cols = {
        "read_id",
        "class_pred",
        "class_probs",
        "barcode_pair"
    }

    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    print("[INFO] Parsing class_probs → P(class 1)...")
    df["p_class1"] = df["class_probs"].apply(extract_p_class1)

    print("[INFO] Applying filters...")
    filtered_df = df[
        (df["class_pred"] == 1) &
        (df["p_class1"] >= PROB_THRESHOLD) &
        (df["barcode_pair"].isin(BARCODE_PAIRS))
    ]

    read_ids = filtered_df["read_id"].unique()

    print(f"[INFO] Writing {len(read_ids)} read IDs to file...")
    os.makedirs(os.path.dirname(OUTPUT_TXT), exist_ok=True)

    with open(OUTPUT_TXT, "w") as f:
        for rid in read_ids:
            f.write(f"{rid}\n")

    print(f"[DONE] Wrote {len(read_ids)} read IDs to:")
    print(f"       {OUTPUT_TXT}")

if __name__ == "__main__":
    main()

