#!/usr/bin/env python3
"""Lipid set enrichment (LipidSig/LION-style) + input files for LION/web.

WHAT THIS DOES
1. Local LSEA: ranks lipids by log2FC vs Untreated and tests whether lipids sharing a
   characteristic (class, category, total carbons, total double bonds, ether, oxidised,
   individual fatty-acyl chains, ...) sit toward the top or bottom of the ranking
   (GSEA-style weighted running sum, permutation p-value, BH FDR). This mirrors the
   "lipid set enrichment" idea in LipidSig / lipidr. It is rank-based, so it works with
   n=1 per condition and is unaffected by a global normalisation offset.
2. LION/web inputs: LION (lipidontology.com) is a web service, so it cannot be run from
   code here. This writes upload-ready files (ranking list, up/down target lists,
   background) with LipidSearch names converted toward LIPID MAPS shorthand, plus a list
   of names LION may not recognise.

WHY NOT LipidSigR HERE: its differential-expression step (deSp_twoGroup, t-test etc.)
needs replicates, and LipidSigR enrichment consumes that result, so it cannot be run on
one sample per condition. Revisit it when you have replicates.

Usage:
    python lipid_enrichment.py lipidomics-data.xlsx [--outdir lipid_results/enrichment]
        [--norm total|median|none] [--min-size 5] [--perms 10000] [--collapse max|none]
"""
import argparse
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lipid_diff_analysis import load, lfc_table, normalise, pseudocount
from lipid_set_analysis import parse, acyl, bh

CATEGORY = {
    **{k: "Glycerophospholipid (GP)" for k in ("PC", "PE", "PS", "PI", "PG", "PA", "CL", "LBPA", "LPC", "LPE", "LPI", "LPG", "LPA", "cPA")},
    **{k: "Glycerolipid (GL)" for k in ("TG", "DG")},
    **{k: "Sphingolipid (SP)" for k in ("Cer", "SM", "SPH", "Hex1Cer", "GM1", "GM2", "ST")},
    **{k: "Sterol (ST)" for k in ("Ch", "ChE", "ZyE")},
    **{k: "Fatty acyl (FA)" for k in ("FA", "AcCa")},
    "Co": "Prenol (PR)",
}
LION_CLASS = {"Hex1Cer": "HexCer", "LBPA": "BMP", "ZyE": "ZyE"}  # LipidSearch -> LIPID MAPS style


def bin_c(c):
    return "<=30" if c <= 30 else "31-34" if c <= 34 else "35-38" if c <= 38 else "39-42" if c <= 42 else ">42"


def bin_db(d):
    return "0-1" if d <= 1 else "2-3" if d <= 3 else "4-5" if d <= 5 else "6-8" if d <= 8 else ">=9"


def features(meta):
    """{lipid ID: set of 'characteristic: value'}"""
    feat = {}
    for i, row in meta.iterrows():
        p = parse(row["Lipid"])
        cls = row["Class"]
        f = {f"class: {cls}", f"category: {CATEGORY.get(cls, 'other')}"}
        ch = [c for c in p["chains"]]
        if ch:
            f.add(f"total carbons: {bin_c(sum(c[1] for c in ch))}")
            f.add(f"total double bonds: {bin_db(sum(c[2] for c in ch))}")
        if p["ether"]:
            f.add("ether-linked: yes")
        if p["oxidised"]:
            f.add("oxidised chain (COOH/CHO): yes")
        if p["molecular"]:
            for c, db in acyl(ch):
                f.add(f"acyl chain: {c}:{db}")
                f.add(f"acyl chain length: {'<16' if c < 16 else '>=24' if c >= 24 else c if c in (16, 18, 20, 22) else 'other'}")
            if any(db >= 4 for _, db in acyl(ch)):
                f.add("PUFA chain (>=4 db): yes")
            if acyl(ch) and all(db == 0 for _, db in acyl(ch)):
                f.add("all acyl chains saturated: yes")
        feat[i] = f
    return feat


def lsea(ranked, sets, perms, rng):
    """Weighted (p=1) GSEA running-sum on `ranked` (Series sorted desc). Returns DataFrame."""
    names = ranked.index.to_numpy()
    r = np.abs(ranked.to_numpy())
    N = len(r)
    pos = {n: k for k, n in enumerate(names)}
    perm_order = np.argsort(rng.random((perms, N)), axis=1)  # shared null draws, first nh columns per set
    rows = []
    for name, members in sets.items():
        idx = np.array(sorted(pos[m] for m in members if m in pos))
        nh = len(idx)
        if nh == 0 or nh >= N:
            continue

        # vectorised: ES for observed and permuted sets
        def es_batch(ixs):  # ixs: (B, nh) positions
            B = ixs.shape[0]
            w = np.zeros((B, N))
            rows_ = np.repeat(np.arange(B), nh)
            w[rows_, ixs.ravel()] = r[ixs].ravel()
            hitc = np.cumsum(w / w.sum(axis=1, keepdims=True), axis=1)
            miss = np.cumsum((w == 0), axis=1) / (N - nh)
            dev = hitc - miss
            mx, mn = dev.max(axis=1), dev.min(axis=1)
            return np.where(mx > -mn, mx, mn)

        obs = es_batch(idx[None, :])[0]
        perm_ix = perm_order[:, :nh]
        null = es_batch(perm_ix)
        same = null[np.sign(null) == np.sign(obs)]
        nes = obs / np.abs(same).mean() if len(same) else np.nan
        p = (np.sum(np.abs(same) >= abs(obs)) + 1) / (len(same) + 1) if len(same) else 1.0
        rows.append(dict(set=name, size=nh, ES=obs, NES=nes, p=p,
                         median_log2FC=np.nan))
    d = pd.DataFrame(rows)
    d["FDR"] = bh(d.p) if len(d) else []
    return d


def lion_name(lipid, cls):
    return re.sub(r"^[^(]+", LION_CLASS.get(cls, cls), lipid, count=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xlsx")
    ap.add_argument("--outdir", default="lipid_results/enrichment")
    ap.add_argument("--sheet", default="RT-Averaged Report")
    ap.add_argument("--norm", default="total", choices=["total", "median", "none"])
    ap.add_argument("--min-size", type=int, default=5)
    ap.add_argument("--perms", type=int, default=10000)
    ap.add_argument("--collapse", default="max", choices=["max", "none"],
                    help="max: keep the highest-area adduct per lipid species (default)")
    ap.add_argument("--lfc", type=float, default=1.0, help="|log2FC| cutoff for LION target lists")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(1)

    mat, meta = load(args.xlsx, args.sheet)
    norm, _ = normalise(mat, args.norm)
    pseudo = pseudocount(norm)
    if args.collapse == "max":
        keep = norm.mean(axis=1).groupby(meta["Lipid"]).idxmax()
        norm, meta = norm.loc[keep], meta.loc[keep]
    feat = features(meta)
    sets = {}
    for i, fs in feat.items():
        for f in fs:
            sets.setdefault(f, []).append(i)
    sets = {k: v for k, v in sets.items() if len(v) >= args.min_size}
    print(f"{len(meta)} lipids, {len(sets)} sets (min size {args.min_size})")

    SHOW = {"Untreated": "Mock", "IFNB": "IFN-β", "IFNG": "IFN-γ"}
    for num, den in (("IFNB", "Untreated"), ("IFNG", "Untreated"), ("IFNG", "IFNB")):
        c = f"{num}_vs_{'Mock' if den == 'Untreated' else den}"      # file tag
        label = f"{SHOW[num]} vs {SHOW[den]}"
        lfc, st = lfc_table(norm, num, den, pseudo)
        both = lfc[st == "both detected"].sort_values(ascending=False)
        # ---------------- LSEA
        sub = {k: [m for m in v if m in both.index] for k, v in sets.items()}
        sub = {k: v for k, v in sub.items() if len(v) >= args.min_size}
        res = lsea(both, sub, args.perms, rng)
        res["median_log2FC"] = [both[sub[s]].median() for s in res.set]
        res = res.sort_values("p")
        res.to_csv(f"{args.outdir}/LSEA_{c}.csv", index=False)
        print(f"\n== {label}: top enriched sets ==")
        print(res.head(12)[["set", "size", "NES", "p", "FDR"]].round(3).to_string(index=False))

        top = res[res.p < 0.05].copy()
        top = top.reindex(top.NES.abs().sort_values(ascending=False).index).head(30).sort_values("NES")
        if len(top):
            col = {"class": "#2471a3", "category": "#16a085", "total carbons": "#8e44ad", "total double bonds": "#d35400",
                   "acyl chain": "#c0392b", "acyl chain length": "#7f8c8d"}
            cc = [col.get(s.split(":")[0], "#555555") for s in top.set]
            fig, ax = plt.subplots(figsize=(8, 0.28 * len(top) + 1.8))
            ax.barh(top.set, top.NES, color=cc, edgecolor=np.where(top.FDR < 0.05, "k", "none"), linewidth=1.2)
            ax.axvline(0, c="k", lw=.6)
            ax.set_xlabel(f"NES (positive = higher in {SHOW[num]})")
            ax.set_title(f"Lipid set enrichment: {label}\n(p<0.05; black outline = FDR<0.05; n=1, exploratory)", fontsize=9)
            ax.tick_params(axis="y", labelsize=7)
            fig.tight_layout(); fig.savefig(f"{args.outdir}/LSEA_{c}.png", dpi=200); plt.close(fig)

        # ---------------- LION/web inputs
        d = pd.DataFrame({"lipid_LION": [lion_name(meta.at[i, "Lipid"], meta.at[i, "Class"]) for i in lfc.index],
                          "lipid_LipidSearch": meta.loc[lfc.index, "Lipid"].values,
                          "log2FC": lfc.values, "detection": st.values}, index=lfc.index)
        d = d[~d.lipid_LipidSearch.str.contains("COOH|CHO")]
        dd = d[d.detection == "both detected"].sort_values("log2FC", ascending=False)
        dd[["lipid_LION", "log2FC"]].to_csv(f"{args.outdir}/LION_{c}_ranking.tsv", sep="\t", index=False, header=False)
        dd[dd.log2FC >= args.lfc].lipid_LION.to_csv(f"{args.outdir}/LION_{c}_target_up.txt", index=False, header=False)
        dd[dd.log2FC <= -args.lfc].lipid_LION.to_csv(f"{args.outdir}/LION_{c}_target_down.txt", index=False, header=False)
        dd.lipid_LION.to_csv(f"{args.outdir}/LION_{c}_background.txt", index=False, header=False)
    meta_ox = meta[meta.Lipid.str.contains("COOH|CHO")].Lipid
    meta_ox.to_csv(f"{args.outdir}/LION_names_excluded_oxidised.txt", index=False, header=False)
    print(f"\nWritten to {args.outdir}/")
    print("LION/web: upload LION_*_ranking.tsv (ranking mode) or the target_up/down + background lists (target-list mode)")


if __name__ == "__main__":
    main()
