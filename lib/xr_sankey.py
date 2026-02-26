#!/usr/bin/env python3
# xr_sankey_pub_v4.py
# ------------------------------------------------------------
# Publication-ready Sankey for Xenoquant pipeline (Plotly version)
# ------------------------------------------------------------
# - Strand-aware demux via BAM flag
# - Clean neutral flows, colored stage nodes
# - Saves vector PDF and CSV summary with losses
# ------------------------------------------------------------

import pandas as pd
import numpy as np
import pysam
from pathlib import Path
import plotly.graph_objects as go
from typing import Set, Dict, Tuple

# ============================================================
# === USER CONFIGURATION =====================================
# ============================================================

TOP_DIR = Path("/top-strand-path")
BOTTOM_DIR = Path("/bottom-strand-path")

SAVE_NAME = "sankey_Xenoquant_pipeline.pdf"
TITLE = "Xenoquant Basecalling — 8L PCR Experiment"

# ============================================================
# === COLOR PALETTE ==========================================
# ============================================================

COLORS = dict(
    pod5   = "#67bfeb",
    bam    = "#6f52a7",
    align  = "#4587b7",
    top    = "#5abe90",
    bottom = "#e94530",
    xna    = "#e9f064",
    dna    = "#cecece",
    loss   = "#b3b3b3",
)

# ============================================================
# === HELPERS ================================================
# ============================================================

def get_pod5_reads(pod5_dir: Path) -> Set[str]:
    from pod5 import Reader
    ids = set()
    for f in pod5_dir.glob("*.pod5"):
        try:
            with Reader(f) as r:
                for read in r:
                    ids.add(read.read_id)
        except Exception:
            pass
    return ids

def get_bam_reads(bam_path: Path) -> Set[str]:
    if not bam_path.exists():
        return set()
    bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
    ids = {r.query_name for r in bam}
    bam.close()
    return ids

def split_strands(bam_path: Path) -> Tuple[Set[str], Set[str]]:
    """Return (top, bottom) read ID sets by strand orientation."""
    top, bottom = set(), set()
    if not bam_path.exists():
        return top, bottom
    bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
    for r in bam:
        if r.is_unmapped:
            continue
        (bottom if r.is_reverse else top).add(r.query_name)
    bam.close()
    return top, bottom

def get_chunks(chunk_path: Path) -> Set[str]:
    if not chunk_path.exists():
        return set()
    npz = np.load(chunk_path, allow_pickle=True)
    for k in ["read_ids", "read_id"]:
        if k in npz:
            return set(npz[k].tolist())
    return set()

def get_demux(csv_path: Path) -> Set[str]:
    if not csv_path.exists():
        return set()
    df = pd.read_csv(csv_path)
    col = [c for c in df.columns if "read" in c.lower()][0]
    return set(df[col])

def get_remora(tsv_path: Path) -> Dict[str, Set[str]]:
    if not tsv_path.exists():
        return {"xna": set(), "dna": set()}
    df = pd.read_csv(tsv_path, sep="\t")
    if "class_pred" not in df.columns or "read_id" not in df.columns:
        return {"xna": set(), "dna": set()}
    return {
        "xna": set(df.loc[df["class_pred"] == 1, "read_id"]),
        "dna": set(df.loc[df["class_pred"] == 0, "read_id"]),
    }

def demux_by_strand(demux_ids: Set[str], bam_path: Path) -> Tuple[Set[str], Set[str]]:
    """Split demuxed reads into top/bottom strands using BAM flag."""
    top, bottom = set(), set()
    if not bam_path.exists() or not demux_ids:
        return top, bottom
    bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
    for r in bam:
        if r.query_name not in demux_ids or r.is_unmapped:
            continue
        (bottom if r.is_reverse else top).add(r.query_name)
    bam.close()
    return top, bottom

# ============================================================
# === COUNT EXTRACTION =======================================
# ============================================================

def collect_counts(top_dir: Path, bottom_dir: Path) -> Dict[str, int]:
    preprocess = top_dir / "preprocess"
    pod5_dir   = preprocess / "pod5"
    bam_dir    = preprocess / "bam"
    bam_path   = bam_dir / "aligned.BAM"

    # --- Stage 1: Raw sets ---
    pod5 = get_pod5_reads(pod5_dir)
    bam  = get_bam_reads(bam_dir / "bc.bam")
    top_aln, bot_aln = split_strands(bam_path)

    # --- Consistency check ---
    total_aln = len(top_aln) + len(bot_aln)
    assert total_aln == len(top_aln) + len(bot_aln), "Strand split mismatch!"

    # --- Stage 2: Chunks ---
    chunks_top    = get_chunks(top_dir / "chunks/basecall_chunks.npz")
    chunks_bottom = get_chunks(bottom_dir / "chunks/basecall_chunks.npz")

    # --- Stage 3: Demux ---
    demux_all = get_demux(top_dir / "demux/all_read_ids.csv")
    chunks_all = chunks_top | chunks_bottom
    demux_all &= chunks_all
    demux_top, demux_bottom = demux_by_strand(demux_all, bam_path)

    # --- Stage 4: Remora outputs ---
    remora_top_raw    = get_remora(top_dir / "remora_outputs/per-read_modifications.tsv")
    remora_bottom_raw = get_remora(bottom_dir / "remora_outputs/per-read_modifications.tsv")

    remora_top = {
        "xna": remora_top_raw["xna"] & demux_top,
        "dna": remora_top_raw["dna"] & demux_top,
    }
    remora_bottom = {
        "xna": remora_bottom_raw["xna"] & demux_bottom,
        "dna": remora_bottom_raw["dna"] & demux_bottom,
    }

    # --- Stage 5: Count summary ---
    counts = dict(
        POD5=len(pod5),
        BAM=len(bam),
        Aligned_Top=len(top_aln),
        Aligned_Bottom=len(bot_aln),
        Chunks_Top=len(chunks_top),
        Chunks_Bottom=len(chunks_bottom),
        Demux_Top=len(demux_top),
        Demux_Bottom=len(demux_bottom),
        Top_XNA=len(remora_top["xna"]),
        Top_DNA=len(remora_top["dna"]),
        Bottom_XNA=len(remora_bottom["xna"]),
        Bottom_DNA=len(remora_bottom["dna"]),
    )

    return counts

# ============================================================
# === PLOT ====================================================
# ============================================================

def make_sankey(counts: Dict[str, int], save_path: Path, title: str, csv_path: Path):
    # --- Labels ---
    labels = [
        "POD5", "BAM", "Aligned BAM",
        "Top strand", "Bottom strand",
        "Top chunks", "Bottom chunks",
        "Top demux", "Bottom demux",
        "Top XNA", "Top DNA", "Bottom XNA", "Bottom DNA",
        "Loss: basecall filtering", "Loss: unaligned",
        "Loss: no chunks", "Loss: undemuxed"
    ]

    colors = [
        COLORS["pod5"], COLORS["bam"], COLORS["align"],
        COLORS["top"], COLORS["bottom"],
        COLORS["top"], COLORS["bottom"],
        COLORS["top"], COLORS["bottom"],
        COLORS["xna"], COLORS["dna"], COLORS["xna"], COLORS["dna"],
        COLORS["loss"], COLORS["loss"], COLORS["loss"], COLORS["loss"],
    ]

    s, t, v, hovertext = [], [], [], []

    def add(src, tgt, val, label):
        if val > 0:
            s.append(src); t.append(tgt); v.append(val)
            hovertext.append(f"{label}<br>{val:,} reads")

    total_aln = counts["Aligned_Top"] + counts["Aligned_Bottom"]

    # Main flows
    add(0, 1, counts["BAM"], "Basecalled reads")
    add(1, 2, total_aln, "Aligned reads")
    add(2, 3, counts["Aligned_Top"], "Top aligned")
    add(2, 4, counts["Aligned_Bottom"], "Bottom aligned")
    add(3, 5, counts["Chunks_Top"], "Top chunks")
    add(4, 6, counts["Chunks_Bottom"], "Bottom chunks")
    add(5, 7, counts["Demux_Top"], "Top demuxed")
    add(6, 8, counts["Demux_Bottom"], "Bottom demuxed")
    add(7, 9, counts["Top_XNA"], "Top XNA")
    add(7, 10, counts["Top_DNA"], "Top DNA")
    add(8, 11, counts["Bottom_XNA"], "Bottom XNA")
    add(8, 12, counts["Bottom_DNA"], "Bottom DNA")

    # --- Loss calculations ---
    basecall_loss = max(counts["POD5"] - counts["BAM"], 0)
    unaligned_loss = max(counts["BAM"] - total_aln, 0)

    nochunk_loss_top = max(counts["Aligned_Top"] - counts["Chunks_Top"], 0)
    nochunk_loss_bottom = max(counts["Aligned_Bottom"] - counts["Chunks_Bottom"], 0)
    nochunk_loss_total = nochunk_loss_top + nochunk_loss_bottom

    undemux_loss_top = max(counts["Chunks_Top"] - counts["Demux_Top"], 0)
    undemux_loss_bottom = max(counts["Chunks_Bottom"] - counts["Demux_Bottom"], 0)
    undemux_loss_total = undemux_loss_top + undemux_loss_bottom

    # Add losses to Sankey
    add(0, 13, basecall_loss, "Basecall filtering loss")
    add(1, 14, unaligned_loss, "Unaligned reads")
    add(3, 15, nochunk_loss_top, "No chunks (top)")
    add(4, 15, nochunk_loss_bottom, "No chunks (bottom)")
    add(5, 16, undemux_loss_top, "Undemuxed top")
    add(6, 16, undemux_loss_bottom, "Undemuxed bottom")



    # --- Build figure ---
    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=30,
            thickness=30,
            line=dict(color="rgba(0,0,0,0)"),
            label=labels,
            color=colors,
        ),
        link=dict(
            source=s,
            target=t,
            value=v,
            color="rgba(0,0,0,0.25)",
            hovertemplate="%{customdata}<extra></extra>",
            customdata=hovertext,
        ),
        textfont=dict(size=14, color="black", family="Arial"),
    ))

    fig.update_layout(
        font=dict(family="Arial", size=14, color="black"),
        title=dict(text=title, font=dict(size=18, family="Arial", color="black")),
        width=1400,
        height=800,
        margin=dict(l=80, r=200, t=80, b=40),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    fig.write_image(str(save_path), scale=3)
    print(f"[Saved] Sankey → {save_path}")

    losses = {
        "Loss_basecall_filtering": basecall_loss,
        "Loss_unaligned": unaligned_loss,
        "Loss_no_chunks_top": nochunk_loss_top,
        "Loss_no_chunks_bottom": nochunk_loss_bottom,
        "Loss_no_chunks_total": nochunk_loss_total,
        "Loss_undemuxed_top": undemux_loss_top,
        "Loss_undemuxed_bottom": undemux_loss_bottom,
        "Loss_undemuxed_total": undemux_loss_total,
    }


    full_stats = {**counts, **losses}
    pd.DataFrame(list(full_stats.items()), columns=["Stage", "Count"]).to_csv(csv_path, index=False)
    print(f"[Saved] Pipeline counts (with losses) → {csv_path}")

# ============================================================
# === MAIN ===================================================
# ============================================================

if __name__ == "__main__":
    counts = collect_counts(TOP_DIR, BOTTOM_DIR)
    print("\n=== PIPELINE COUNTS ===")
    for k, v in counts.items():
        print(f"{k:20s}: {v:,}")

    sankey_dir = TOP_DIR / "sankey"
    sankey_dir.mkdir(exist_ok=True)

    save_path = sankey_dir / SAVE_NAME
    csv_path  = sankey_dir / "pipeline_counts.csv"

    make_sankey(counts, save_path, TITLE, csv_path)

