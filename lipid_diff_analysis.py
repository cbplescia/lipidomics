#!/usr/bin/env python3
"""Differential lipid abundance for the pilot: IFNB vs Untreated, IFNG vs Untreated,
and how similar the two IFN responses are.

Samples (one per condition, n=1): s5 = Untreated, s6 = IFN-beta, s7 = IFN-gamma.

IMPORTANT: with a single sample per group there is no variance estimate, so this
script reports effect sizes (log2 fold changes) only -- no p-values. Treat the
output as hypothesis generating and use it to design the replicated experiment.

Usage:
    python lipid_diff_analysis.py lipidomics-data.xlsx [--outdir results]
        [--sheet "RT-Averaged Report"] [--norm total|median|none] [--lfc 1.0]
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

SAMPLES = {"s5": "Untreated", "s6": "IFNB", "s7": "IFNG"}
HEADER_ROW = 3  # 0-based row of the column header in the sheet


def load(path, sheet):
    raw = pd.read_excel(path, sheet_name=sheet, header=HEADER_ROW)
    area = {f"Raw Area — {s}": c for s, c in SAMPLES.items()}
    missing = [c for c in area if c not in raw.columns]
    if missing:
        raise SystemExit(f"Columns not found: {missing}\nFound: {list(raw.columns)}")

    df = raw.copy()
    # footer / note rows have no Class or mass
    df = df[df["Class"].notna() & df["Calc. Mass (Da)"].notna()].copy()

    name = df["Lipid Species"].fillna(df["Lipid Molecular Group"])
    df["Lipid"] = name
    df["ID"] = name.astype(str) + " [" + df["Adduct"].astype(str) + "]"
    df["ID"] = df["ID"].where(~df["ID"].duplicated(keep=False),
                              df["ID"] + " RT" + df["RT (min)"].round(2).astype(str))

    # internal standards: names like 'PC(15:0_18:1)+D7:(s)' and 'Ch-D7'
    is_std = df["Lipid"].astype(str).str.contains(r"\+D\d+:\(s\)|-D\d+$", regex=True)
    df["is_IS"] = is_std
    print(f"Loaded {len(df)} rows; removing {is_std.sum()} internal standards")

    df = df.rename(columns=area)
    df = df[~df["is_IS"]].set_index("ID")
    keep = list(area.values())
    mat = df[keep].apply(pd.to_numeric, errors="coerce")
    mat.columns = list(SAMPLES.values())
    mat = mat.mask(mat <= 0)  # 0 / blank = not detected
    meta = df[["Lipid", "Class", "Adduct", "RT (min)", "Grade"]]
    return mat, meta


def normalise(mat, how):
    if how == "none":
        return mat, pd.Series(1.0, index=mat.columns)
    if how == "total":
        size = mat.sum()
    else:  # median of log ratios vs. the mean profile, on lipids seen in all samples
        full = mat.dropna()
        ref = np.log2(full).mean(axis=1)
        size = 2 ** (np.log2(full).sub(ref, axis=0).median())
    size = size / size.mean()
    return mat / size, size


def lfc_table(mat, num, den, pseudo):
    a, b = mat[num], mat[den]
    lfc = np.log2((a.fillna(0) + pseudo) / (b.fillna(0) + pseudo))
    status = pd.Series("both detected", index=mat.index)
    status[a.isna() & b.notna()] = f"only in {den}"
    status[a.notna() & b.isna()] = f"only in {num}"
    return lfc, status


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xlsx")
    ap.add_argument("--outdir", default="lipid_results")
    ap.add_argument("--sheet", default="RT-Averaged Report")
    ap.add_argument("--norm", default="total", choices=["total", "median", "none"],
                    help="between-sample normalisation (default: total signal)")
    ap.add_argument("--lfc", type=float, default=1.0,
                    help="|log2FC| threshold for calling a change (default 1 = 2-fold)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mat, meta = load(args.xlsx, args.sheet)
    norm, size = normalise(mat, args.norm)
    print(f"{len(mat)} lipids; size factors ({args.norm}):\n{size.round(3).to_string()}")
    # pseudocount = half the smallest detected value, so "not detected" stays finite
    pseudo = norm.min().min() / 2

    out = meta.copy()
    out[[f"{c}_area" for c in norm.columns]] = norm.values
    for cond in ["IFNB", "IFNG"]:
        out[f"log2FC_{cond}_vs_Untreated"], out[f"detect_{cond}_vs_Untreated"] = \
            lfc_table(norm, cond, "Untreated", pseudo)
    out["log2FC_IFNG_vs_IFNB"], out["detect_IFNG_vs_IFNB"] = \
        lfc_table(norm, "IFNG", "IFNB", pseudo)

    x, y = out["log2FC_IFNB_vs_Untreated"], out["log2FC_IFNG_vs_Untreated"]
    t = args.lfc

    def call(v):
        return np.where(v >= t, "up", np.where(v <= -t, "down", "ns"))

    cb, cg = call(x), call(y)
    pattern = np.select(
        [(cb == "ns") & (cg == "ns"),
         (cb != "ns") & (cg == "ns"),
         (cb == "ns") & (cg != "ns"),
         cb == cg],
        ["unchanged", "IFNB-specific", "IFNG-specific", "shared (same direction)"],
        default="opposite direction")
    out["IFNB_call"], out["IFNG_call"], out["response_pattern"] = cb, cg, pattern
    out["delta_log2FC (IFNG-IFNB)"] = y - x
    out = out.sort_values("delta_log2FC (IFNG-IFNB)", key=abs, ascending=False)
    out.to_csv(f"{args.outdir}/all_lipids_differential.csv")
    for cond in ["IFNB", "IFNG"]:
        col = f"log2FC_{cond}_vs_Untreated"
        out.sort_values(col, key=abs, ascending=False).to_csv(
            f"{args.outdir}/{cond}_vs_Untreated.csv")

    # ---- similarity of the two responses ----
    ok = out.dropna(subset=["log2FC_IFNB_vs_Untreated", "log2FC_IFNG_vs_Untreated"])
    xb, yg = ok["log2FC_IFNB_vs_Untreated"], ok["log2FC_IFNG_vs_Untreated"]
    pr, ps = stats.pearsonr(xb, yg), stats.spearmanr(xb, yg)
    chg = ok[(ok.IFNB_call != "ns") | (ok.IFNG_call != "ns")]
    lines = [
        f"Lipids analysed: {len(ok)} (normalisation: {args.norm}; threshold |log2FC|>={t})",
        f"IFNB vs Untreated: {(ok.IFNB_call=='up').sum()} up, {(ok.IFNB_call=='down').sum()} down",
        f"IFNG vs Untreated: {(ok.IFNG_call=='up').sum()} up, {(ok.IFNG_call=='down').sum()} down",
        "",
        "Response pattern counts:",
        ok.response_pattern.value_counts().to_string(),
        "",
        f"Correlation of log2FC(IFNB) vs log2FC(IFNG), all lipids: "
        f"Pearson r={pr[0]:.3f}, Spearman rho={ps[0]:.3f}",
    ]
    if len(chg) > 3:
        pc = stats.pearsonr(chg.log2FC_IFNB_vs_Untreated, chg.log2FC_IFNG_vs_Untreated)
        sign = (np.sign(chg.log2FC_IFNB_vs_Untreated) ==
                np.sign(chg.log2FC_IFNG_vs_Untreated)).mean()
        lines.append(f"Among {len(chg)} lipids changed in either: Pearson r={pc[0]:.3f}, "
                     f"same direction for {sign:.0%}")
    lines.append("\nNOTE: n=1 per condition -> no p-values; fold changes include "
                 "technical + biological noise.")
    summary = "\n".join(lines)
    print("\n" + summary)
    open(f"{args.outdir}/summary.txt", "w").write(summary + "\n")

    # per-class view
    cls = ok.groupby("Class").agg(
        n=("Lipid", "size"),
        median_log2FC_IFNB=("log2FC_IFNB_vs_Untreated", "median"),
        median_log2FC_IFNG=("log2FC_IFNG_vs_Untreated", "median"))
    cls.sort_values("n", ascending=False).to_csv(f"{args.outdir}/class_summary.csv")

    # ---- plots ----
    colors = {"unchanged": "#bbbbbb", "IFNB-specific": "#1f77b4",
              "IFNG-specific": "#d62728", "shared (same direction)": "#2ca02c",
              "opposite direction": "#ff7f0e"}
    fig, ax = plt.subplots(figsize=(7, 6.5))
    for p, c in colors.items():
        s = ok[ok.response_pattern == p]
        ax.scatter(s.log2FC_IFNB_vs_Untreated, s.log2FC_IFNG_vs_Untreated, s=14, c=c,
                   alpha=.75, label=f"{p} ({len(s)})", edgecolor="none")
    lim = np.nanmax(np.abs(ok[["log2FC_IFNB_vs_Untreated", "log2FC_IFNG_vs_Untreated"]].values)) * 1.05
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=.7)
    ax.axhline(0, c="k", lw=.4); ax.axvline(0, c="k", lw=.4)
    for v in (-t, t):
        ax.axhline(v, c="grey", lw=.4, ls=":"); ax.axvline(v, c="grey", lw=.4, ls=":")
    top = ok.reindex(ok["delta_log2FC (IFNG-IFNB)"].abs().sort_values(ascending=False).index[:8])
    top = pd.concat([top, ok.reindex((xb.abs() + yg.abs()).sort_values(ascending=False).index[:8])]).drop_duplicates()
    for i, r in top.iterrows():
        ax.annotate(r.Lipid, (r.log2FC_IFNB_vs_Untreated, r.log2FC_IFNG_vs_Untreated),
                    fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2FC IFNB vs Untreated"); ax.set_ylabel("log2FC IFNG vs Untreated")
    ax.set_title(f"IFNB vs IFNG response (Pearson r={pr[0]:.2f}, Spearman={ps[0]:.2f})")
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    fig.tight_layout(); fig.savefig(f"{args.outdir}/IFNB_vs_IFNG_log2FC_scatter.png", dpi=200)
    plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, cond in zip(axs, ["IFNB", "IFNG"]):
        lf = ok[f"log2FC_{cond}_vs_Untreated"]
        avg = np.log10(ok[["Untreated", cond]].mean(axis=1) if False else
                       (ok[f"{cond}_area"] + ok["Untreated_area"]).div(2).fillna(
                           ok[f"{cond}_area"].fillna(ok["Untreated_area"])))
        sig = lf.abs() >= t
        ax.scatter(avg[~sig], lf[~sig], s=10, c="#bbbbbb", edgecolor="none")
        ax.scatter(avg[sig], lf[sig], s=14, c="#d62728", edgecolor="none")
        for i in lf.abs().sort_values(ascending=False).index[:8]:
            ax.annotate(ok.at[i, "Lipid"], (avg[i], lf[i]), fontsize=6,
                        xytext=(3, 3), textcoords="offset points")
        ax.axhline(0, c="k", lw=.4)
        for v in (-t, t):
            ax.axhline(v, c="grey", lw=.4, ls=":")
        ax.set_xlabel("log10 mean normalised area"); ax.set_title(f"{cond} vs Untreated")
    axs[0].set_ylabel("log2 fold change")
    fig.tight_layout(); fig.savefig(f"{args.outdir}/MA_plots.png", dpi=200); plt.close(fig)

    # heatmap of the most changed lipids (log2 vs Untreated)
    n = 40
    idx = (xb.abs() + yg.abs()).sort_values(ascending=False).index[:n]
    h = ok.loc[idx, ["log2FC_IFNB_vs_Untreated", "log2FC_IFNG_vs_Untreated"]]
    h.columns = ["IFNB", "IFNG"]
    fig, ax = plt.subplots(figsize=(4.5, 0.22 * n + 1.5))
    m = np.nanmax(np.abs(h.values))
    im = ax.imshow(h.values, cmap="RdBu_r", vmin=-m, vmax=m, aspect="auto")
    ax.set_xticks([0, 1]); ax.set_xticklabels(h.columns)
    ax.set_yticks(range(n)); ax.set_yticklabels(ok.loc[idx, "Lipid"], fontsize=6)
    fig.colorbar(im, label="log2FC vs Untreated"); ax.set_title(f"Top {n} changed lipids")
    fig.tight_layout(); fig.savefig(f"{args.outdir}/top_changed_heatmap.png", dpi=200)
    plt.close(fig)

    # class-level bar plot
    c2 = cls[cls.n >= 3].sort_values("n", ascending=False)
    fig, ax = plt.subplots(figsize=(9, 4))
    w = .4; xs = np.arange(len(c2))
    ax.bar(xs - w/2, c2.median_log2FC_IFNB, w, label="IFNB", color="#1f77b4")
    ax.bar(xs + w/2, c2.median_log2FC_IFNG, w, label="IFNG", color="#d62728")
    ax.set_xticks(xs); ax.set_xticklabels([f"{k}\n(n={v})" for k, v in c2.n.items()], fontsize=7)
    ax.axhline(0, c="k", lw=.5); ax.set_ylabel("median log2FC vs Untreated"); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{args.outdir}/class_median_log2FC.png", dpi=200); plt.close(fig)

    print(f"\nResults written to {args.outdir}/")


if __name__ == "__main__":
    main()
