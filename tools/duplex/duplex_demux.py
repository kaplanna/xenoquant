import pandas as pd

# === File Paths (edit these as needed) ===
demux_csv_path = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/BA-Basecall/demux/all_read_ids.csv"
duplex_txt_path = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/duplex_read_IDs.txt"
output_csv_path = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/duplex_read_matches.csv"

# === Load demux file ===
demux_df = pd.read_csv(demux_csv_path)

# === Parse duplex read pairs (dx:1 only) ===
duplex_pairs = []
with open(duplex_txt_path, "r") as f:
    in_duplex_section = False
    for line in f:
        line = line.strip()
        if line.startswith("dx:1"):
            in_duplex_section = True
            continue
        elif line.startswith("dx:0"):
            break  # Stop when reaching unpaired reads
        if in_duplex_section and ";" in line:
            read_id_1, read_id_2 = line.split(";")
            duplex_pairs.append((read_id_1, read_id_2))

# === Convert to DataFrame ===
duplex_df = pd.DataFrame(duplex_pairs, columns=["read_id_1", "read_id_2"])

# === Merge both read_ids with demux info ===
merged_1 = duplex_df.merge(demux_df, left_on="read_id_1", right_on="read_id", how="left").rename(
    columns={"sample_id": "sample_id_1", "barcode_pair": "barcode_pair_1"})
merged_2 = merged_1.merge(demux_df, left_on="read_id_2", right_on="read_id", how="left").rename(
    columns={"sample_id": "sample_id_2", "barcode_pair": "barcode_pair_2"})

# Drop redundant columns
merged_2.drop(columns=["read_id_x", "read_id_y"], inplace=True)

# === Add flag for same barcode pair ===
merged_2["same_barcode_pair"] = merged_2["barcode_pair_1"] == merged_2["barcode_pair_2"]

# === Save merged results ===
merged_2.to_csv(output_csv_path, index=False)

# === Summary Stats ===
same_count = merged_2["same_barcode_pair"].sum()
total_pairs = len(merged_2)
different_count = total_pairs - same_count

print("Duplex Pair Summary:")
print(f"Total duplex pairs:        {total_pairs}")
print(f"Same barcode pair:         {same_count}")
print(f"Different barcode pair:    {different_count}")

