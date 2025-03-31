import pandas as pd

# === File paths ===
duplex_csv = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/duplex_read_matches.csv"
ba_file = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/BA-Basecall/demux/demux_per-read_modifications.tsv"
st_file = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/ST-Basecall/demux/demux_per-read_modifications.tsv"
output_file = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/duplex_predictions.csv"
summary_file = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/duplex_summary.txt"
valid_output_file = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/duplex_valid_only.csv"

# === Load data ===
duplex_df = pd.read_csv(duplex_csv)
ba_df = pd.read_csv(ba_file, sep="\t", usecols=["read_id", "class_pred"])
st_df = pd.read_csv(st_file, sep="\t", usecols=["read_id", "class_pred"])

# === Set index for lookup ===
ba_dict = ba_df.set_index("read_id")["class_pred"].to_dict()
st_dict = st_df.set_index("read_id")["class_pred"].to_dict()

# === Assign predictions ===
duplex_df["BA_pred"] = duplex_df["read_id_1"].map(ba_dict).combine_first(
                        duplex_df["read_id_2"].map(ba_dict))
duplex_df["ST_pred"] = duplex_df["read_id_1"].map(st_dict).combine_first(
                        duplex_df["read_id_2"].map(st_dict))

# === Determine source of each read ===
duplex_df["read1_source"] = duplex_df["read_id_1"].map(lambda r: "BA" if r in ba_dict else "ST" if r in st_dict else "None")
duplex_df["read2_source"] = duplex_df["read_id_2"].map(lambda r: "BA" if r in ba_dict else "ST" if r in st_dict else "None")

def classify_duplex_type(r1, r2):
    combo = f"{r1}-{r2}"
    if combo in ["BA-ST", "ST-BA", "BA-BA", "ST-ST", "None-None"]:
        return combo
    else:
        return "Ambiguous"

duplex_df["duplex_type"] = duplex_df.apply(
    lambda row: classify_duplex_type(row["read1_source"], row["read2_source"]), axis=1
)

# === Count duplex types ===
duplex_type_counts = duplex_df["duplex_type"].value_counts()

print("\n📦 Duplex Type Breakdown:")
for dtype, count in duplex_type_counts.items():
    print(f"{dtype}: {count}")

# === Filter to valid duplex types only ===
valid_duplex_df = duplex_df[duplex_df["duplex_type"].isin(["BA-ST", "ST-BA"])].copy()
# === Further filter: keep only duplexes from the same barcode pair ===
valid_duplex_df = valid_duplex_df[valid_duplex_df["same_barcode_pair"] == True].copy()

valid_duplex_df.to_csv(valid_output_file, index=False)

# === Recalculate stats ONLY from valid filtered data ===
ba_count = valid_duplex_df["BA_pred"].notna().sum()
st_count = valid_duplex_df["ST_pred"].notna().sum()
both_count = valid_duplex_df[(valid_duplex_df["BA_pred"].notna()) & (valid_duplex_df["ST_pred"].notna())].shape[0]
neither_count = ((valid_duplex_df["BA_pred"].isna()) & (valid_duplex_df["ST_pred"].isna())).sum()

both_1 = ((valid_duplex_df["BA_pred"] == 1) & (valid_duplex_df["ST_pred"] == 1)).sum()
both_0 = ((valid_duplex_df["BA_pred"] == 0) & (valid_duplex_df["ST_pred"] == 0)).sum()
ba1_st0 = ((valid_duplex_df["BA_pred"] == 1) & (valid_duplex_df["ST_pred"] == 0)).sum()
st1_ba0 = ((valid_duplex_df["BA_pred"] == 0) & (valid_duplex_df["ST_pred"] == 1)).sum()

summary = f"""
📊 Duplex Prediction Summary (Valid Only: BA-ST/ST-BA + Same Barcode)

Read Source Breakdown:
-----------------------
Reads with BA_pred:        {ba_count}
Reads with ST_pred:        {st_count}
Reads in BOTH (expected):  {both_count}
Reads in neither:          {neither_count}

Prediction Combinations:
------------------------
Both 1 (BA=1, ST=1):        {both_1}
Both 0 (BA=0, ST=0):        {both_0}
BA=1, ST=0:                 {ba1_st0}
ST=1, BA=0:                 {st1_ba0}

Final Output:
------------------------
Valid duplex pairs kept (BA-ST + ST-BA, same barcode): {len(valid_duplex_df)}
"""
# === Breakdown by barcode pair ===
barcode_summary_lines = []
barcode_pairs = valid_duplex_df["barcode_pair_1"].unique()

barcode_summary_lines.append("\n🔬 Breakdown by Barcode Pair:\n")
for barcode in sorted(barcode_pairs):
    sub_df = valid_duplex_df[valid_duplex_df["barcode_pair_1"] == barcode]

    total = len(sub_df)
    b1_s1 = ((sub_df["BA_pred"] == 1) & (sub_df["ST_pred"] == 1)).sum()
    b0_s0 = ((sub_df["BA_pred"] == 0) & (sub_df["ST_pred"] == 0)).sum()
    b1_s0 = ((sub_df["BA_pred"] == 1) & (sub_df["ST_pred"] == 0)).sum()
    b0_s1 = ((sub_df["BA_pred"] == 0) & (sub_df["ST_pred"] == 1)).sum()

    barcode_summary_lines.append(f"{barcode} — Total: {total} | Both 1: {b1_s1} | Both 0: {b0_s0} | BA=1/ST=0: {b1_s0} | ST=1/BA=0: {b0_s1}")

# Print to terminal
print("\n".join(barcode_summary_lines))

# Also append to summary file
with open(summary_file, "a") as f:
    f.write("\n".join(barcode_summary_lines))

print(summary)

with open(summary_file, "w") as f:
    f.write("📦 Duplex Type Breakdown:\n")
    for dtype, count in duplex_type_counts.items():
        f.write(f"{dtype}: {count}\n")
    f.write(summary)

print(f"✅ Annotated duplex file saved to {output_file}")
print(f"✅ Valid filtered duplex file saved to {valid_output_file}")
print(f"📝 Summary saved to {summary_file}")

