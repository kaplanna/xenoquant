import pandas as pd

# === File Paths (edit these as needed) ===
ba_file = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/BA-Basecall/demux/demux_per-read_modifications.tsv"
st_file = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/ST-Basecall/demux/demux_per-read_modifications.tsv"

import pandas as pd

def count_shared_read_ids(ba_file, st_file):
    # Load TSV files
    df1 = pd.read_csv(ba_file, sep='\t')
    df2 = pd.read_csv(st_file, sep='\t')

    # Extract read_ids
    read_ids1 = set(df1['read_id'])
    read_ids2 = set(df2['read_id'])

    # Find shared read_ids
    shared_ids = read_ids1.intersection(read_ids2)

    # Output results
    print(f"BA file: {len(read_ids1)} read_ids")
    print(f"ST file: {len(read_ids2)} read_ids")
    print(f"Shared read_ids: {len(shared_ids)}")

    return shared_ids

shared_read_ids = count_shared_read_ids(ba_file, st_file)

