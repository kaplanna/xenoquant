import pandas as pd
import os

# === Input file paths ===
GC_file = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/PZn/250328_P2-P3_Basecall/PG-Basecall/demux/demux_per-read_modifications.tsv'
AT_file = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/PZn/250604_P2-P3_Basecall/PA-Basecall/demux/demux_per-read_modifications.tsv'

# === Output directory ===
output_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/PZn/250604_P2-P3_Basecall/PGA_model_comparison'
os.makedirs(output_dir, exist_ok=True)
summary_output = os.path.join(output_dir, 'model_comparison_by_barcode_pair.csv')
detailed_output = os.path.join(output_dir, 'per_read_classification.csv')

# === Load input files ===
df_GC = pd.read_csv(GC_file, sep='\t')
df_AT = pd.read_csv(AT_file, sep='\t')

def split_probs(probs_str):
    try:
        p0, p1 = map(float, probs_str.split(','))
        return pd.Series({'class_0_probs': p0, 'class_1_probs': p1})
    except Exception:
        return pd.Series({'class_0_probs': None, 'class_1_probs': None})

# Apply to both GC and AT
df_GC[['class_0_probs', 'class_1_probs']] = df_GC['class_probs'].apply(split_probs)
df_AT[['class_0_probs', 'class_1_probs']] = df_AT['class_probs'].apply(split_probs)


# === Merge on read_id ===
merged_df = pd.merge(df_GC, df_AT, on='read_id', suffixes=('_GC', '_AT'))


# === Classification logic ===
def classify(row):
    at_pred = row['class_pred_AT']
    gc_pred = row['class_pred_GC']
    
    if at_pred == 1 and gc_pred == 1:
        return 'XNA'
    elif at_pred == 1 and gc_pred == 0:
        return 'GC'
    elif at_pred == 0 and gc_pred == 1:
        return 'AT'
    else:
        return 'AT' if row['class_0_probs_AT'] > row['class_0_probs_GC'] else 'GC'

merged_df['final_call'] = merged_df.apply(classify, axis=1)

# === Select and save detailed output for inspection ===
detailed_cols = [
    'read_id',
    'barcode_pair_GC',
    'class_pred_AT', 'class_0_probs_AT',
    'class_pred_GC', 'class_0_probs_GC',
    'final_call'
]
detailed_df = merged_df[detailed_cols]
detailed_df.to_csv(detailed_output, index=False)
print(f"Saved per-read classification to: {detailed_output}")

# === Summarize by barcode_pair ===
summary = detailed_df.groupby(['barcode_pair_GC', 'final_call']).size().unstack(fill_value=0)
summary['total'] = summary.sum(axis=1)
summary['XNA_fraction'] = summary.get('XNA', 0) / summary['total']
summary['AT_fraction'] = summary.get('AT', 0) / summary['total']
summary['GC_fraction'] = summary.get('GC', 0) / summary['total']

# === Save summary output ===
summary.to_csv(summary_output)
print(f"Saved barcode summary to: {summary_output}")
print(summary)

