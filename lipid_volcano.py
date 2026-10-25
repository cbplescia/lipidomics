#!/usr/bin/env python3
"""Volcano-style plots: IFNB vs Untreated and IFNG vs Untreated.

With one sample per condition there is no variance, so a p-value (the usual
volcano y-axis) cannot be computed and none is invented here. Instead the
y-axis is mean abundance of the two samples (log10 normalised area): a fold
change at high abundance is better measured than one near the detection limit.
Hollow triangles = detected in only one of the two samples (value is capped by
the pseudocount floor).

Usage:
    python lipid_volcano.py lipidomics-data.xlsx [--outdir lipid_results/volcano]
        [--norm total|median|none] [--lfc 1.0] [--label 15]
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lipid_diff_analysis import load, lfc_table, normalise, pseudocount

COND = {"IFNB": "IFN-β", "IFNG": "IFN-γ"}
UP, DOWN, NS = "#c0392b", "#2471a3", "#bbbbbb"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xlsx")
    ap.add_argument("--outdir", default="lipid_results/volcano")
    ap.add_argument("--sheet", default="RT-Averaged Report")
    ap.add_argument("--norm", default="total", choices=["total", "median", "none"])
    ap.add_argument("--lfc", type=float, default=1.0)
    ap.add_argument("--label", type=int, default=15, help="labels per plot")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mat, meta = load(args.xlsx, args.sheet)
    norm, _ = normalise(mat, args.norm)
    pseudo = pseudocount(norm)
    t = args.lfc

    res = {}
    for c in COND:
        lfc, st = lfc_table(norm, c, "Untreated", pseudo)
        d = pd.DataFrame({"lfc": lfc, "status": st}).join(meta)
        # mean of the two samples, not-detected counted as 0
        d["y"] = np.log10(norm[[c, "Untreated"]].fillna(0).mean(axis=1))
        res[c] = d
    lim = max(res[c].lfc.abs().max() for c in COND) * 1.08
    ymax = max(res[c].y.max() for c in COND) * 1.03
    ymin = min(res[c].y.min() for c in COND) * 0.97

    for c, label in COND.items():
        d = res[c]
        fig, ax = plt.subplots(figsize=(7.5, 6.5))
        one = d.status != "both detected"
        col = np.where(d.lfc >= t, UP, np.where(d.lfc <= -t, DOWN, NS))
        ax.scatter(d.lfc[~one], d.y[~one], c=col[~one.values], s=16, alpha=.8,
                   edgecolor="none")
        ax.scatter(d.lfc[one], d.y[one], c="none", edgecolor=col[one.values], marker="^",
                   s=26, linewidth=.9)
        for v in (-t, t):
            ax.axvline(v, c="grey", lw=.5, ls=":")
        ax.axvline(0, c="k", lw=.5)
        # label strongest changes weighted toward well-measured lipids
        score = d.lfc.abs() * (d.y - ymin)
        for i in score[d.lfc.abs() >= t].sort_values(ascending=False).index[:args.label]:
            ax.annotate(d.at[i, "Lipid"], (d.at[i, "lfc"], d.at[i, "y"]), fontsize=6.5,
                        xytext=(3, 3), textcoords="offset points")
        ax.set_xlim(-lim, lim); ax.set_ylim(ymin, ymax)
        ax.set_xlabel(f"log2 fold change vs Untreated ({args.norm} normalised)")
        ax.set_ylabel("log10 mean normalised area (Untreated, treated)")
        nu, nd = (d.lfc >= t).sum(), (d.lfc <= -t).sum()
        ax.set_title(f"{label} vs Untreated: {nu} up, {nd} down of {len(d)}  (n=1, no p-values)",
                     fontsize=10)
        fig.tight_layout()
        fig.savefig(f"{args.outdir}/volcano_{c}_vs_Untreated.png", dpi=200)
        plt.close(fig)
    print(f"Plots written to {args.outdir}/")


if __name__ == "__main__":
    main()
