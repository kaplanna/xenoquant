# xr_model_paths.py

# Maps forward strand mod base to its reverse complement
mod_base_pairs = {
    "B": "S",
    "P": "Z",
    "Ds": "Px"
}

# All models (keyed by context on B-strand, as written in filename)
model_paths = {
    "GBA": {
        "BA": "/home/marchandlab/github/kaplanna/xemora/models/240930_NTC_Models/GBA-BA-model_best.pt",
        "ST": "/home/marchandlab/github/kaplanna/xemora/models/240930_NTC_Models/GBA-ST-model_best.pt"
    },
    "GBT": {
        "BA": "/home/marchandlab/github/kaplanna/xemora/models/240930_NTC_Models/GBT-BA-model_best.pt",
        "ST": "/home/marchandlab/github/kaplanna/xemora/models/240930_NTC_Models/GBT-ST-model_best.pt"
    },
    "GDsA": {
        "DsA": "/home/marchandlab/github/kaplanna/xemora/models/241212_DsPx_Models/Ds-Px/G-Ds-A/Ds-A-model_best.pt",
        "PxT": "/home/marchandlab/github/kaplanna/xemora/models/241212_DsPx_Models/Ds-Px/G-Ds-A/Px-T-model_best.pt"
    },
    "TPC": {
        "PG": "/home/marchandlab/github/kaplanna/xemora/models/250321_PZ_Models/PG-model_best.pt",
        "ZC": "/home/marchandlab/github/kaplanna/xemora/models/250321_PZ_Models/ZC-model_best.pt"
    }
}

