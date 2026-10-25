#!/usr/bin/env python3
"""Per-lipid-group log2FC bar plots, separately for IFNB and IFNG vs Untreated.

Usage:
    python lipid_group_plots.py lipidomics-data.xlsx [--outdir lipid_results/group_plots]
        [--norm total|median|none]

Groups (by LipidSearch Class):
  PIPs           PI, LPI   (no phosphorylated PIPs were detected in this dataset)
  PS             PS
  Ceramides      Cer, Hex1Cer
  Sphingolipids  SM, SPH, GM1, GM2, ST
  Cholesterols   Ch, ChE, ZyE (no eicosanoids/oxylipins were detected)
Bars are log2 fold change vs Untreated (normalised area). Hatched, outlined bars
mean the lipid was not detected in one of the two samples (value uses a pseudocount).
n=1 per condition: effect sizes only, no statistics.
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lipid_diff_analysis import load, lfc_table, normalise, pseudocount

GROUPS = {
    "PI_LPI": ["PI", "LPI"],
    "PS": ["PS"],
    "Ceramides": ["Cer", "Hex1Cer"],
    "Sphingolipids": ["SM", "SPH", "GM1", "GM2", "ST"],
    "Cholesterols": ["Ch", "ChE", "ZyE"],
}
COND = {"IFNB": "IFN-β", "IFNG": "IFN-γ"}
UP, DOWN = "#c0392b", "#2471a3"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xlsx")
    ap.add_argument("--outdir", default="lipid_results/group_plots")
    ap.add_argument("--sheet", default="RT-Averaged Report")
    ap.add_argument("--norm", default="total", choices=["total", "median", "none"])
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mat, meta = load(args.xlsx, args.sheet)
    norm, _ = normalise(mat, args.norm)
    pseudo = pseudocount(norm)
    res = {}
    for c in COND:
        lfc, st = lfc_table(norm, c, "Untreated", pseudo)
        res[c] = pd.DataFrame({"lfc": lfc, "status": st}).join(meta)

    # same x range for IFNB and IFNG so the two are visually comparable
    lim = max(res[c].lfc.abs().max() for c in COND) * 1.05

    for g, classes in GROUPS.items():
        for c, label in COND.items():
            d = res[c][res[c].Class.isin(classes)].sort_values("lfc")
            if d.empty:
                continue
            fig, ax = plt.subplots(figsize=(7, 1.2 + 0.2 * len(d)))
            y = np.arange(len(d))
            nd = d.status != "both detected"
            ax.barh(y, d.lfc, color=np.where(d.lfc >= 0, UP, DOWN),
                    edgecolor=np.where(nd, "k", "none"), linewidth=.8)
            for yi in y[nd.values]:
                ax.patches[yi].set_hatch("////")
            ax.set_yticks(y)
            ax.set_yticklabels(d.Lipid, fontsize=7)
            ax.axvline(0, c="k", lw=.6)
            for v in (-1, 1):
                ax.axvline(v, c="grey", lw=.5, ls=":")
            ax.set_xlim(-lim, lim)
            ax.set_xlabel(f"log2 fold change vs Untreated ({args.norm} normalised)")
            ax.set_title(f"{g.replace('_', ' ')} — {label} vs Untreated  "
                         f"({(d.lfc>=1).sum()} up, {(d.lfc<=-1).sum()} down of {len(d)})",
                         fontsize=9)
            ax.set_title(ax.get_title() + "\nhatched = not detected in one sample; dotted = 2-fold",
                         fontsize=9)
            fig.tight_layout()
            fig.savefig(f"{args.outdir}/{g}_{c}_vs_Untreated.png", dpi=200)
            plt.close(fig)
            d.to_csv(f"{args.outdir}/{g}_{c}_vs_Untreated.csv")
    print(f"Plots written to {args.outdir}/")


if __name__ == "__main__":
    main()
