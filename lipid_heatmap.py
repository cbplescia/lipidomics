#!/usr/bin/env python3
"""Heatmap of every lipid species: Mock (Untreated) vs IFN-γ (optionally IFN-β too).

Rows   : all lipid species, grouped by class (classes ordered by total abundance), and
         sorted within each class by log2FC (IFN-γ vs Mock, highest at top).
Left   : class colour strip.
Panel 1: log10 normalised area in each sample (grey = not detected).
Panel 2: log2 fold change vs Mock (diverging colours, capped at --cap).
Strip  : black tick = detected in only one of the two samples (log2FC is then set by the
         not-detected floor, so treat as presence/absence).

n=1 per condition: descriptive only. Normalisation default = total signal.

Usage:
    python lipid_heatmap.py lipidomics-data.xlsx [--outdir lipid_results/heatmap]
        [--norm total|median|none] [--include-ifnb] [--classes PC PE TG ...]
        [--cap 4] [--min-class 1] [--row-height 0.055]
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

from lipid_diff_analysis import load, lfc_table, normalise, pseudocount

SHOW = {"Untreated": "Mock", "IFNB": "IFN-β", "IFNG": "IFN-γ"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xlsx")
    ap.add_argument("--outdir", default="lipid_results/heatmap")
    ap.add_argument("--sheet", default="RT-Averaged Report")
    ap.add_argument("--norm", default="total", choices=["total", "median", "none"])
    ap.add_argument("--include-ifnb", action="store_true", help="add IFN-β columns")
    ap.add_argument("--classes", nargs="+", help="only these LipidSearch classes")
    ap.add_argument("--cap", type=float, default=4.0, help="log2FC colour limit (default ±4)")
    ap.add_argument("--row-height", type=float, default=0.055, help="inches per species (lower = shorter figure)")
    ap.add_argument("--min-class", type=int, default=1, help="drop classes with fewer species")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mat, meta = load(args.xlsx, args.sheet)
    norm, _ = normalise(mat, args.norm)
    pseudo = pseudocount(norm)
    conds = ["IFNG"] + (["IFNB"] if args.include_ifnb else [])
    samples = ["Untreated"] + (["IFNB"] if args.include_ifnb else []) + ["IFNG"]

    df = meta[["Lipid", "Class"]].copy()
    for s in samples:
        df[f"area_{s}"] = norm[s]
    for c in conds:
        df[f"lfc_{c}"], df[f"status_{c}"] = lfc_table(norm, c, "Untreated", pseudo)
    if args.classes:
        df = df[df.Class.isin(args.classes)]
    df = df[df.groupby("Class").Class.transform("size") >= args.min_class]

    # order: classes by total abundance, species by IFN-γ log2FC (high at top)
    corder = df.groupby("Class")["area_Untreated"].sum().add(df.groupby("Class")["area_IFNG"].sum(), fill_value=0)
    corder = corder.sort_values(ascending=False).index.tolist()
    df["_c"] = df.Class.map({c: i for i, c in enumerate(corder)})
    df = df.sort_values(["_c", "lfc_IFNG"], ascending=[True, False])
    df.drop(columns="_c").to_csv(f"{args.outdir}/heatmap_order_and_values.csv")

    n = len(df)
    A = np.log10(df[[f"area_{s}" for s in samples]].to_numpy(float))
    L = df[[f"lfc_{c}" for c in conds]].to_numpy(float)
    one = np.column_stack([(df[f"status_{c}"] != "both detected").to_numpy() for c in conds])
    A = np.ma.masked_invalid(A)

    # class colours
    pal = plt.get_cmap("tab20")
    ccol = {c: pal(i % 20) for i, c in enumerate(corder)}
    strip = np.array([ccol[c][:3] for c in df.Class])[:, None, :]

    wA, wL = len(samples), len(conds)
    fig = plt.figure(figsize=(3.4 + 0.9 * wA + 0.9 * wL, max(6, args.row_height * n + 2)))
    gs = fig.add_gridspec(1, 5, width_ratios=[0.25, wA, wL, 0.25 * wL, 0.02], wspace=0.06,
                          left=0.17, right=0.9, top=1 - 1.4 / fig.get_figheight(), bottom=0.9 / fig.get_figheight())
    axs = [fig.add_subplot(gs[0, i]) for i in range(4)]
    axs[0].imshow(strip, aspect="auto", interpolation="nearest")
    cm1 = plt.get_cmap("viridis").copy(); cm1.set_bad("#eeeeee")
    vmin, vmax = np.nanpercentile(A.compressed(), [2, 99.5])
    im1 = axs[1].imshow(A, aspect="auto", cmap=cm1, vmin=vmin, vmax=vmax, interpolation="nearest")
    cm2 = plt.get_cmap("RdBu_r")
    im2 = axs[2].imshow(np.clip(L, -args.cap, args.cap), aspect="auto", cmap=cm2, vmin=-args.cap, vmax=args.cap, interpolation="nearest")
    axs[3].imshow(one.astype(float), aspect="auto", cmap=ListedColormap(["white", "black"]), interpolation="nearest")

    # class labels / separators
    cls_arr = df.Class.to_numpy()
    start = 0
    for c in corder:
        k = int((cls_arr == c).sum())
        if k == 0:
            continue
        for ax in axs[:4]:
            ax.axhline(start - 0.5, c="k", lw=0.4)
        axs[0].text(-0.7, start + k / 2 - 0.5, f"{c} ({k})", ha="right", va="center", fontsize=7 if k > 6 else 6)
        start += k
    for ax in axs:
        ax.set_yticks([]); ax.set_xticks([])
    axs[1].set_xticks(range(wA)); axs[1].set_xticklabels([SHOW[s] for s in samples], fontsize=8)
    axs[2].set_xticks(range(wL)); axs[2].set_xticklabels([f"{SHOW[c]}\nvs Mock" for c in conds], fontsize=8)
    axs[3].set_xticks(range(wL)); axs[3].set_xticklabels(["1-sample\nonly"] * wL, fontsize=6)
    for ax in axs[1:]:
        ax.xaxis.tick_top()
    axs[1].set_title("log10 normalised area", fontsize=8, pad=22)
    axs[2].set_title("log2 fold change", fontsize=8, pad=22)
    cax1 = fig.add_axes([0.17, 0.35 / fig.get_figheight(), 0.18, 0.12 / fig.get_figheight() * 2])
    fig.colorbar(im1, cax=cax1, orientation="horizontal").set_label("log10 area", fontsize=7)
    cax2 = fig.add_axes([0.5, 0.35 / fig.get_figheight(), 0.18, 0.12 / fig.get_figheight() * 2])
    fig.colorbar(im2, cax=cax2, orientation="horizontal", extend="both").set_label(f"log2FC (±{args.cap:g} cap)", fontsize=7)
    for c in (cax1, cax2):
        c.tick_params(labelsize=6)
    fig.suptitle(f"All lipid species: Mock vs IFN-γ{' (and IFN-β)' if args.include_ifnb else ''}  "
                 f"({n} species, {args.norm}-normalised; n=1, descriptive)", fontsize=9, y=0.995)
    for ext in ("png", "pdf"):
        fig.savefig(f"{args.outdir}/lipid_species_heatmap.{ext}", dpi=200)
    plt.close(fig)
    print(f"{n} species plotted -> {args.outdir}/lipid_species_heatmap.png (+ .pdf, + order CSV)")


if __name__ == "__main__":
    main()
