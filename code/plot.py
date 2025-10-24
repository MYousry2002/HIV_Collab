#!/usr/bin/env python3
"""
cluster_heatmap.py

Read a count matrix (TSV) whose first column is 'sequence_name' and
plot a clustered heatmap with both row and column dendrograms.

- Handles wide label sets cleanly (bigger canvas, rotated x labels)
- Optional log1p transform
- Optional z-score scaling by rows or columns
- Uses correlation distance + average linkage by default (good for TF profiles)
- Filters zero-variance rows/columns (after transform/scale) to avoid clustering errors
- Saves PNG and PDF

Usage:
  python plot.py \
      --input ../results/fimo_hiv1_clean/counts.tsv \
      --output ../results/fimo_hiv1_clean/plots/heatmap.png \
      --pdf \
      --log1p \
      --zscore rows \
      --metric correlation \
      --method average \
      --cmap viridis
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--input", required=True, help="TSV with first column 'sequence_name'.")
    p.add_argument("--output", required=True, help="Output image path (PNG).")
    p.add_argument("--pdf", action="store_true", help="Also write a PDF with same basename.")
    p.add_argument("--log1p", action="store_true", help="Apply log1p transform before scaling.")
    p.add_argument("--zscore", choices=["none", "rows", "cols"], default="none",
                   help="Z-score by rows or columns (after log1p if requested).")
    p.add_argument("--metric", default="correlation",
                   help="Distance metric for clustering (e.g., correlation, euclidean, cosine).")
    p.add_argument("--method", default="average",
                   help="Linkage method (e.g., average, ward, complete, single).")
    p.add_argument("--cmap", default="viridis",
                   help="Matplotlib/Seaborn colormap (e.g., viridis, magma, vlag, coolwarm).")
    p.add_argument("--figwidth", type=float, default=20, help="Figure width (inches).")
    p.add_argument("--figheight", type=float, default=22, help="Figure height (inches).")
    p.add_argument("--dpi", type=int, default=300, help="Output DPI for PNG.")
    p.add_argument("--xtick_rot", type=int, default=75, help="Rotation for column labels.")
    p.add_argument("--ytick_size", type=float, default=7.5, help="Font size for row labels.")
    p.add_argument("--dendro_ratio_rows", type=float, default=0.15, help="Row dendrogram height ratio.")
    p.add_argument("--dendro_ratio_cols", type=float, default=0.08, help="Col dendrogram width ratio.")
    p.add_argument("--cbar", action="store_true", help="Show colorbar.")
    p.add_argument("--clip_percentile", type=float, default=None,
                   help="If set (e.g., 99), clip heatmap values to +/- this percentile of |values| "
                        "(good for outlier-heavy counts).")
    p.add_argument("--exclude-col", type=str, default=None, help="Column name to exclude from the heatmap (e.g., 'ZNF').")
    return p.parse_args()


def zscore_df(df: pd.DataFrame, axis: str) -> pd.DataFrame:
    if axis == "rows":
        mu = df.mean(axis=1)
        sd = df.std(axis=1).replace(0, np.nan)
        return (df.sub(mu, axis=0)).div(sd, axis=0)
    elif axis == "cols":
        mu = df.mean(axis=0)
        sd = df.std(axis=0).replace(0, np.nan)
        return (df - mu) / sd
    else:
        return df


def filter_zero_variance(df: pd.DataFrame) -> pd.DataFrame:
    # Drop rows/cols with all-NaN or zero variance
    # (after transform/scale these can appear).
    row_var = df.var(axis=1, numeric_only=True)
    col_var = df.var(axis=0, numeric_only=True)
    keep_rows = row_var.replace(0, np.nan).notna()
    keep_cols = col_var.replace(0, np.nan).notna()
    filtered = df.loc[keep_rows, keep_cols]
    return filtered


def compute_linkage(df: pd.DataFrame, metric: str, method: str):
    # pdist on values; returns condensed distance vector
    # Handle the single-row/col edge cases gracefully.
    if df.shape[0] > 1:
        row_dist = pdist(df.values, metric=metric)
        row_link = linkage(row_dist, method=method)
    else:
        row_link = None

    if df.shape[1] > 1:
        col_dist = pdist(df.values.T, metric=metric)
        col_link = linkage(col_dist, method=method)
    else:
        col_link = None

    return row_link, col_link


def main():
    args = parse_args()

    # ---- Load ----
    df = pd.read_csv(args.input, sep="\t", header=0)
    if df.columns[0] != "sequence_name":
        # Be forgiving: treat the first column as index anyway
        df = df.rename(columns={df.columns[0]: "sequence_name"})
    df = df.set_index("sequence_name")
    if args.exclude_col and args.exclude_col in df.columns:
        df = df.drop(columns=[args.exclude_col])
        print(f"[info] Excluded column: {args.exclude_col}")

    # Ensure numeric (coerce errors to NaN)
    df = df.apply(pd.to_numeric, errors="coerce")

    # ---- Transform / Scale ----
    mat = df.copy()
    if args.log1p:
        mat = np.log1p(mat)

    if args.zscore in {"rows", "cols"}:
        mat = zscore_df(mat, "rows" if args.zscore == "rows" else "cols")

    # ---- Filter zero-variance / all-NaN ----
    mat = mat.dropna(how="all", axis=0).dropna(how="all", axis=1)
    mat = filter_zero_variance(mat)

    if mat.shape[0] == 0 or mat.shape[1] == 0:
        raise ValueError("After filtering, the matrix is empty. Check inputs or disable aggressive filtering.")

    # ---- Optional clipping for outliers (nice for raw counts) ----
    vmin = vmax = None
    if args.clip_percentile is not None:
        p = float(args.clip_percentile)
        if args.zscore == "none":
            # symmetric clipping around median can be odd for counts; just clip high end
            vmax = np.nanpercentile(mat.values, p)
            vmin = np.nanpercentile(mat.values, 100 - p)
        else:
            # z-scored: symmetric is fine
            a = np.nanpercentile(np.abs(mat.values), p)
            vmin, vmax = -a, a

    # ---- Compute dendrograms explicitly to avoid version mismatches ----
    row_link, col_link = compute_linkage(mat, metric=args.metric, method=args.method)

    # ---- Plot ----
    sns.set(context="notebook")
    # Note: Using clustermap so dendrograms are aligned and sized predictably
    g = sns.clustermap(
        mat,
        row_linkage=row_link,
        col_linkage=col_link,
        cmap=args.cmap,
        vmin=vmin,
        vmax=vmax,
        cbar=args.cbar,
        figsize=(args.figwidth, args.figheight),
        dendrogram_ratio=(args.dendro_ratio_rows, args.dendro_ratio_cols),
        xticklabels=True,
        yticklabels=True
    )

    # Label aesthetics
    # Rotate x tick labels and shrink y tick labels to reduce overlap
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=args.xtick_rot, ha="right")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=args.ytick_size)

    # Put colorbar on the right (if present) and tighten layout
    if args.cbar and g.cax is not None:
        g.cax.set_title("value", fontsize=9)

    plt.tight_layout()

    # ---- Save ----
    out_png = Path(args.output)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    g.figure.savefig(out_png, dpi=args.dpi, bbox_inches="tight")

    if args.pdf:
        out_pdf = out_png.with_suffix(".pdf")
        g.figure.savefig(out_pdf, bbox_inches="tight")

    print(f"[ok] saved: {out_png}")
    if args.pdf:
        print(f"[ok] saved: {out_pdf}")


if __name__ == "__main__":
    main()