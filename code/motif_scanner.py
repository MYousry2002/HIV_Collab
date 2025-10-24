#!/usr/bin/env python3
"""
Scan sequences with FIMO using a MEME-format motif file and produce:
  (1) all_hits.tsv  : raw FIMO hits (fimo.tsv passthrough with minor cleanup)
  (2) counts.tsv    : sequence x motif count matrix (drop all-zero motifs)

Requires:
  - MEME Suite (fimo) on PATH
  - Python 3.8+ with pandas

Example:
  python motif_scanner.py \
    --meme ../data/nonredundunt_tf_motifs.meme \
    --fasta ../data/HIV-1_LTR.fasta \
    --outdir ../results/fimo_hiv1 \
    --thresh 1e-4

  python motif_scanner.py \
    --meme ../data/nonredundunt_tf_motifs.meme \
    --fasta ../data/HIV-2_LTR.fasta \
    --outdir ../results/fimo_hiv2 \
    --thresh 1e-4

  # with metadata to map clusters -> TF names
  python motif_scanner.py \
    --meme ../data/nonredundunt_tf_motifs.meme \
    --fasta ../data/HIV-1_LTR.fasta \
    --meta  ../data/motif_metadata.tsv \
    --outdir ../results/fimo_hiv1 \
    --thresh 1e-4

  # manual label cleaning (no metadata):
  python motif_scanner.py \
    --meme ../data/nonredundunt_tf_motifs.meme \
    --fasta ../data/HIV-1_LTR.fasta \
    --outdir ../results/fimo_hiv1_clean \
    --thresh 1e-4 \
    --clean-labels \
    --family-map ../data/TF_family.csv
    
"""
from __future__ import annotations
import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import itertools
import collections
import csv
from typing import List

def load_metadata(meta_path: Path) -> pd.DataFrame:
    """
    Load metadata TSV expected to contain at least columns: 'cluster' and 'tf_name'.
    Additional columns (e.g., motif_id, source_id, family_name) are ignored here.
    """
    dfm = pd.read_csv(meta_path, sep="\t", dtype=str)
    required = {"cluster", "tf_name"}
    missing = required - set(dfm.columns)
    if missing:
        sys.exit(f"[!] Metadata is missing required columns: {missing}")
    # Normalize cluster identifiers (e.g., ensure they look like AC####)
    dfm["cluster"] = dfm["cluster"].str.strip()
    dfm["tf_name"] = dfm["tf_name"].str.strip()
    return dfm

def build_cluster_to_tf_map(dfm: pd.DataFrame) -> dict:
    """
    Build a mapping: cluster (e.g., AC0001) -> representative TF name string.
    If multiple TF names are present for a cluster, join unique names (up to 3) with '|',
    appending '…' if there are more.
    """
    mapping = {}
    for clus, sub in dfm.groupby("cluster"):
        names = [x for x in sub["tf_name"].dropna().astype(str).str.strip().unique() if x]
        if not names:
            continue
        names_sorted = sorted(names)
        if len(names_sorted) <= 3:
            label = "|".join(names_sorted)
        else:
            label = "|".join(names_sorted[:3]) + "|…"
        mapping[clus] = label
    return mapping

def extract_cluster_from_motif_id(motif_id: str) -> str | None:
    """Extract leading cluster like 'AC0001' from FIMO motif_id such as 'AC0001:GATA/PROP:GATA'."""
    if not isinstance(motif_id, str):
        return None
    m = re.match(r"^(AC\d{4})", motif_id)
    return m.group(1) if m else None

def clean_motif_label(motif_id: str) -> str | None:
    """Manually clean a motif_id like 'AC0389:FOXN:Forkhead' -> 'FOXN';
    'AC0390:CREB/JUND:bZIP' -> 'CREB/JUND'. If no colon after the AC tag, return the residual token.
    """
    if not isinstance(motif_id, str):
        return None
    s = motif_id.strip()
    # Drop leading cluster prefix if present
    s = re.sub(r"^AC\d{4}:", "", s)
    # Keep text up to the next ':' (if any)
    if ":" in s:
        s = s.split(":", 1)[0]
    # Remove any bracketed annotations lingering in name
    s = re.sub(r"\[.*?\]|\(.*?\)", "", s).strip()
    return s or None

def aggregate_counts_by_clean_label(mat_by_motif: pd.DataFrame) -> pd.DataFrame:
    """Aggregate a motif_id-count matrix into cleaned labels using clean_motif_label()."""
    if mat_by_motif.empty:
        return mat_by_motif
    # Map each motif_id column to a cleaned label
    col_map = {}
    for col in mat_by_motif.columns:
        label = clean_motif_label(col)
        if label:
            col_map[col] = label
    if not col_map:
        return pd.DataFrame(index=mat_by_motif.index)
    # Group columns by cleaned label and sum
    groups: dict[str, list[str]] = {}
    for col, lbl in col_map.items():
        groups.setdefault(lbl, []).append(col)
    out = pd.DataFrame(index=mat_by_motif.index)
    for lbl, cols in groups.items():
        out[lbl] = mat_by_motif[cols].sum(axis=1).astype(int)
    # Drop all-zero columns and sort
    if out.shape[1] > 0:
        out = out.loc[:, (out.sum(axis=0) > 0)]
        out = out.sort_index(axis=0).sort_index(axis=1)
    return out

def normalize_composite_label(name: str) -> str:
    """Alphabetize slash-separated composite labels: 'KLF/SP' <-> 'SP/KLF'."""
    if not isinstance(name, str):
        return name
    name = name.strip()
    if "/" in name:
        parts = [p.strip() for p in name.split("/") if p.strip()]
        parts.sort()
        return "/".join(parts)
    return name

def load_family_map(csv_path: Path) -> dict:
    """Read TF_family mapping CSV with columns: TF_family,motif_clusters.
    Returns dict mapping cleaned+normalized motif label -> TF_family.
    """
    dfm = pd.read_csv(csv_path)
    required = {"TF_family", "motif_clusters"}
    missing = required - set(dfm.columns)
    if missing:
        sys.exit(f"[!] Family map CSV missing required columns: {missing}")
    label_to_family: dict[str, str] = {}
    for _, row in dfm.iterrows():
        fam = str(row["TF_family"]).strip()
        clusters = str(row["motif_clusters"]).split(",")
        for raw in clusters:
            lab = raw.strip()
            if not lab:
                continue
            # Clean like our motif labels and normalize composites
            lab = normalize_composite_label(lab)
            label_to_family[lab] = fam
    return label_to_family

def aggregate_counts_by_family(mat_by_motif: pd.DataFrame, label_to_family: dict) -> pd.DataFrame:
    """Aggregate motif counts into TF families using the provided label->family map.
    Uses cleaned+normalized labels as the key to map columns; unmapped labels are kept as-is.
    """
    if mat_by_motif.empty:
        return mat_by_motif
    # First derive cleaned labels from motif_id columns
    col_map_clean: dict[str, str] = {}
    for col in mat_by_motif.columns:
        cl = clean_motif_label(col)
        if cl is None:
            continue
        cl = normalize_composite_label(cl)
        col_map_clean[col] = cl
    # Now map cleaned labels to families (fallback to the cleaned label if not in map)
    col_to_family: dict[str, str] = {}
    for col, cl in col_map_clean.items():
        fam = label_to_family.get(cl, cl)  # keep original cleaned label if not mapped
        col_to_family[col] = fam
    # Group columns by family and sum
    groups: dict[str, list[str]] = {}
    for col, fam in col_to_family.items():
        groups.setdefault(fam, []).append(col)
    out = pd.DataFrame(index=mat_by_motif.index)
    for fam, cols in groups.items():
        out[fam] = mat_by_motif[cols].sum(axis=1).astype(int)
    # Drop all-zero columns and sort
    if out.shape[1] > 0:
        out = out.loc[:, (out.sum(axis=0) > 0)]
        out = out.sort_index(axis=0).sort_index(axis=1)
    return out


# Deduplicate overlapping motif hits **across all motifs/labels** within each sequence (and strand if present)
def deduplicate_overlaps(df: pd.DataFrame, overlap_frac: float = 0.3) -> pd.DataFrame:
    """Deduplicate overlapping hits **across all motifs/labels** within each sequence.
    If two hits overlap by more than `overlap_frac` of the shorter length, keep only the
    higher-scoring hit. Deduplication is done per (sequence_name[, strand]).
    - Assumes numeric columns: start, stop, score.
    - If 'strand' exists, dedup per strand separately.
    """
    if df.empty:
        return df
    need = {"sequence_name", "start", "stop", "score"}
    missing = need - set(df.columns)
    if missing:
        return df  # required columns not present; return unchanged

    per_strand = "strand" in df.columns
    out_rows: List[pd.Series] = []

    group_cols = ["sequence_name"] + (["strand"] if per_strand else [])
    for _, sub in df.groupby(group_cols, sort=False):
        if len(sub) == 1:
            out_rows.append(sub.iloc[0])
            continue
        sub_sorted = sub.sort_values(by=["score"], ascending=False).reset_index(drop=True)
        kept_intervals: list[tuple[int, int]] = []
        kept_rows: list[pd.Series] = []
        for _, row in sub_sorted.iterrows():
            s1 = int(row["start"]); e1 = int(row["stop"])  # ensure ordered
            if e1 < s1:
                s1, e1 = e1, s1
            l1 = e1 - s1 + 1
            overlaps = False
            for (ks, ke) in kept_intervals:
                ov = min(e1, ke) - max(s1, ks) + 1
                if ov > 0:
                    l2 = ke - ks + 1
                    shorter = l1 if l1 < l2 else l2
                    if ov > overlap_frac * shorter:
                        overlaps = True
                        break
            if not overlaps:
                kept_intervals.append((s1, e1))
                kept_rows.append(row)
        if kept_rows:
            out_rows.extend(kept_rows)
    if not out_rows:
        return df.iloc[0:0].copy()
    return pd.DataFrame(out_rows).reset_index(drop=True)

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meme", required=True, help="MEME-format motif file (archetypes).")
    ap.add_argument("--fasta", required=True, help="FASTA of sequences to scan.")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--thresh", type=float, default=1e-4, help="FIMO p-value threshold (default: 1e-4).")
    ap.add_argument("--bg", choices=["uniform"], default="uniform",
                    help="Background model (currently only 'uniform' is applied).")
    ap.add_argument("--fimo-bin", default="fimo", help="Path to fimo executable (default: fimo in PATH).")
    ap.add_argument("--meta", help="TSV metadata with columns: cluster, tf_name (and others). If provided, outputs will include cluster and tf_name; counts will be grouped by tf_name.")
    ap.add_argument("--clean-labels", action="store_true",
                    help="Manually clean motif_id column labels (drop AC####: prefix; keep token before the next colon). If set, counts will be grouped by these cleaned labels.")
    ap.add_argument("--family-map", dest="family_map", help="CSV mapping TF_family,motif_clusters. If provided, counts will be aggregated by TF_family and written to counts_by_family.tsv and counts.tsv.")
    return ap.parse_args()

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def clean_fasta(in_fa: Path, out_fa: Path):
    """
    Uppercase; replace anything not A/C/G/T with N (so '-' or other IUPAC letters don't crash FIMO).
    """
    valid = set("ACGT")
    with open(in_fa, "r") as inf, open(out_fa, "w") as outf:
        seq_id = None
        buf = []
        def flush():
            if seq_id is not None:
                outf.write(f">{seq_id}\n")
                # wrap at 80 chars for readability
                s = "".join(buf)
                for i in range(0, len(s), 80):
                    outf.write(s[i:i+80] + "\n")
        for line in inf:
            if line.startswith(">"):
                flush()
                seq_id = line[1:].strip()
                buf = []
            else:
                s = line.strip().upper()
                # remove alignment dashes entirely so motifs aren't broken
                s = s.replace("-", "")
                # convert any remaining non-ACGT to N
                s = "".join(ch if ch in valid else "N" for ch in s)
                buf.append(s)
        flush()

def run_fimo(fimo_bin: str, meme: Path, cleaned_fa: Path, oc_dir: Path, pthresh: float):
    """
    Run fimo. Produces oc_dir/fimo.tsv among other files.
    """
    # If the output directory exists from a previous run, remove to avoid "non-empty" error
    if oc_dir.exists():
        shutil.rmtree(oc_dir)
    cmd = [
        fimo_bin,
        "--thresh", str(pthresh),
        "--verbosity", "1",
        # add --bgfile if you later want non-uniform; default is uniform which matches our MEME header
        "--oc", str(oc_dir),
        str(meme),
        str(cleaned_fa),
    ]
    print("[fimo] Running:", " ".join(cmd), file=sys.stderr)
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"[!] FIMO failed with exit code {e.returncode}")

def read_fimo_tsv(fimo_dir: Path) -> pd.DataFrame:
    fimo_tsv = fimo_dir / "fimo.tsv"
    if not fimo_tsv.exists():
        sys.exit(f"[!] Expected FIMO output not found: {fimo_tsv}")
    # FIMO uses tab-delimited with a leading comment line starting with '#'
    df = pd.read_csv(fimo_tsv, sep="\t", comment="#", dtype=str)
    # Cast numeric cols where appropriate if present
    for col in ["start", "stop", "score", "p-value", "q-value"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df

def write_all_hits(df: pd.DataFrame, out_path: Path):
    # Keep original FIMO columns; just write out
    df.to_csv(out_path, sep="\t", index=False)

def make_count_matrix(df: pd.DataFrame, label_col: str = "motif_id") -> pd.DataFrame:
    """
    Build counts: number of hits per (sequence_name, <label_col>), where label_col is either
    'motif_id' (default) or 'tf_name' (when metadata provided). Drops columns with all-zero sums.
    """
    required = {"sequence_name", label_col}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        sys.exit(f"[!] FIMO table missing required columns for counting: {missing}")

    grp = df.groupby(["sequence_name", label_col], dropna=False).size().rename("count").reset_index()
    mat = grp.pivot(index="sequence_name", columns=label_col, values="count").fillna(0).astype(int)
    # Drop all-zero columns
    if mat.shape[1] > 0:
        mat = mat.loc[:, (mat.sum(axis=0) > 0)]
    # Sort rows/cols for stability
    mat = mat.sort_index(axis=0)
    if mat.shape[1] > 0:
        mat = mat.sort_index(axis=1)
    return mat

def main():
    args = parse_args()
    outdir = Path(args.outdir)
    ensure_dir(outdir)

    meme = Path(args.meme)
    fasta = Path(args.fasta)
    if not meme.exists():
        sys.exit(f"[!] MEME file not found: {meme}")
    if not fasta.exists():
        sys.exit(f"[!] FASTA not found: {fasta}")

    # Optional metadata mapping (cluster -> tf_name)
    meta_path = Path(args.meta) if getattr(args, "meta", None) else None
    df_meta = None
    cluster_to_tf = {}
    if meta_path is not None:
        if not meta_path.exists():
            sys.exit(f"[!] Metadata TSV not found: {meta_path}")
        print(f"[i] Loading metadata -> {meta_path}", file=sys.stderr)
        df_meta = load_metadata(meta_path)
        cluster_to_tf = build_cluster_to_tf_map(df_meta)

    # 1) Clean FASTA (handle '-' and other non-ACGT as N)
    cleaned_fa = outdir / "cleaned.fa"
    print(f"[i] Cleaning FASTA -> {cleaned_fa}", file=sys.stderr)
    clean_fasta(fasta, cleaned_fa)

    # 2) Run FIMO
    fimo_dir = outdir / "fimo_out"
    run_fimo(args.fimo_bin, meme, cleaned_fa, fimo_dir, args.thresh)

    # 3) Read FIMO hits
    df = read_fimo_tsv(fimo_dir)

    # If metadata provided, add cluster + tf_name columns derived from motif_id
    if cluster_to_tf:
        df["cluster"] = df["motif_id"].apply(extract_cluster_from_motif_id)
        df["tf_name"] = df["cluster"].map(cluster_to_tf)

    # Clean label always (used for family mapping / manual cleaning)
    df["clean_label"] = df["motif_id"].apply(clean_motif_label)
    df["clean_label"] = df["clean_label"].apply(normalize_composite_label)

    # Load family map early (if provided) to deduplicate at the FAMILY level
    label_to_family = None
    dedup_label_col = "motif_id"  # default
    if getattr(args, "family_map", None):
        fam_path = Path(args.family_map)
        if not fam_path.exists():
            sys.exit(f"[!] Family map CSV not found: {fam_path}")
        print(f"[i] Loading TF family map -> {fam_path}", file=sys.stderr)
        label_to_family = load_family_map(fam_path)
        # Map cleaned labels to family; fallback to cleaned label if not mapped
        df["family_label"] = df["clean_label"].map(lambda x: label_to_family.get(x, x))
        dedup_label_col = "family_label"
    elif args.clean_labels:
        dedup_label_col = "clean_label"
    elif cluster_to_tf:
        dedup_label_col = "tf_name"

    # Deduplicate overlapping hits **across motifs** within each sequence (and per strand if present)
    before = len(df)
    df = deduplicate_overlaps(df, overlap_frac=0.5)
    after = len(df)
    print(f"[i] Global dedup of overlapping hits (>50% of shorter motif) per sequence/strand: {before} -> {after} hits", file=sys.stderr)

    # 4) Write all hits (post-dedup; includes optional columns)
    all_hits = outdir / "all_hits.tsv"
    print(f"[i] Writing all hits -> {all_hits}", file=sys.stderr)
    write_all_hits(df, all_hits)

    # 5) Count matrix
    if len(df) == 0:
        print("[i] No FIMO hits found at the given threshold; counts files will be empty after column drop.", file=sys.stderr)
        (outdir / "counts_by_motif.tsv").write_text("")
        (outdir / "counts_by_tf.tsv").write_text("")
        (outdir / "counts.tsv").write_text("")
        print("[done] all_hits.tsv and empty counts written.", file=sys.stderr)
        return

    # Always create motif-level counts
    mat_by_motif = make_count_matrix(df, label_col="motif_id")
    counts_by_motif_path = outdir / "counts_by_motif.tsv"
    print(f"[i] Writing count matrix by motif_id -> {counts_by_motif_path}", file=sys.stderr)
    mat_by_motif.to_csv(counts_by_motif_path, sep="\t")

    wrote_primary = False
    # (A) If a TF family mapping CSV is provided, aggregate to families (primary output)
    if getattr(args, "family_map", None):
        if label_to_family is None:
            sys.exit("[!] Internal error: family map not loaded earlier.")
        mat_by_family = aggregate_counts_by_family(mat_by_motif, label_to_family)
        counts_by_family_path = outdir / "counts_by_family.tsv"
        print(f"[i] Writing count matrix by TF family -> {counts_by_family_path}", file=sys.stderr)
        mat_by_family.to_csv(counts_by_family_path, sep="\t")
        # counts.tsv mirrors the family-aggregated matrix when provided
        mat_by_family.to_csv(outdir / "counts.tsv", sep="\t")
        wrote_primary = True

    # (B) Else, manual cleaned labels requested -> aggregate by cleaned label
    if (not wrote_primary) and args.clean_labels:
        mat_by_clean = aggregate_counts_by_clean_label(mat_by_motif)
        counts_by_tf_path = outdir / "counts_by_tf.tsv"
        print(f"[i] Writing count matrix by cleaned label -> {counts_by_tf_path}", file=sys.stderr)
        mat_by_clean.to_csv(counts_by_tf_path, sep="\t")
        mat_by_clean.to_csv(outdir / "counts.tsv", sep="\t")
        wrote_primary = True

    # (C) Else, if metadata provided, aggregate by TF name via cluster mapping
    if (not wrote_primary) and cluster_to_tf:
        mat_by_tf = aggregate_counts_by_tf(mat_by_motif, cluster_to_tf)
        counts_by_tf_path = outdir / "counts_by_tf.tsv"
        print(f"[i] Writing count matrix by TF name -> {counts_by_tf_path}", file=sys.stderr)
        mat_by_tf.to_csv(counts_by_tf_path, sep="\t")
        mat_by_tf.to_csv(outdir / "counts.tsv", sep="\t")
        wrote_primary = True

    # (D) Fallback: mirror counts_by_motif
    if not wrote_primary:
        mat_by_motif.to_csv(outdir / "counts.tsv", sep="\t")

    print("[done] all_hits.tsv and counts.tsv written.", file=sys.stderr)

if __name__ == "__main__":
    main()