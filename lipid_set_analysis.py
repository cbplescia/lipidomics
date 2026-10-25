#!/usr/bin/env python3
"""Lipid "pathway / GO-style" set analysis for IFNB and IFNG vs Untreated.

Each lipid is assigned to one or more functional sets (membrane role, organelle,
signalling, storage, chemical features such as ether-linked / oxidised / PUFA-
containing). Set-level shifts are then summarised per contrast and drawn as a
GO-style dot plot and a per-set distribution plot.

Statistics: for each set, a Mann-Whitney test compares the log2FC of the set's
lipids against all other lipids (BH-adjusted across sets). "Relative log2FC" is
the set median minus the median of all lipids, so a global normalisation offset
(e.g. less total lipid in one sample) cancels out. Lipids are treated as
independent and there is one sample per condition, so this is EXPLORATORY: it
ranks which functional groups move together, it does not prove significance.
Lipids detected in only one of the two samples are excluded from set statistics.

Usage:
    python lipid_set_analysis.py lipidomics-data.xlsx [--outdir lipid_results/sets]
        [--norm total|median|none] [--min-size 4]
Edit SETS below to add or change groupings.
"""
import argparse
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from lipid_diff_analysis import load, lfc_table, normalise, pseudocount

# ---------------------------------------------------------------- parsing
CHAIN = re.compile(r"^(?P<pre>[A-Za-z]*-?)(?P<c>\d+):(?P<db>\d+)(?P<mod>[A-Za-z]*)$")


def parse(name):
    """Return dict(cls, chains[(pre,C,db,mod)], molecular, ether, oxidised)."""
    m = re.match(r"^([^(]+)\((.*)\)$", str(name))
    if not m:
        return dict(cls=str(name), chains=[], molecular=False, ether=False, oxidised=False)
    cls, inner = m.groups()
    toks = inner.split("_")
    chains = []
    for t in toks:
        c = CHAIN.match(t)
        if c:
            chains.append((c["pre"], int(c["c"]), int(c["db"]), c["mod"]))
    return dict(
        cls=cls, chains=chains, molecular=len(toks) > 1,
        ether=any(p in ("O-", "P-") for p, *_ in chains),
        oxidised=any(mod in ("COOH", "CHO", "OH", "OOH") for *_, mod in chains))


def acyl(chains):  # fatty-acyl chains only (not sphingoid bases, not oxidised)
    return [(c, db) for pre, c, db, mod in chains if pre in ("", "O-", "P-") and not mod]


def has_pufa(p, n=4):
    return p["molecular"] and any(db >= n for _, db in acyl(p["chains"]))


def total_db(p):
    return sum(db for *_, db, _m in p["chains"])


# ---------------------------------------------------------------- set definitions
# (set name, category, biological role, membership function(lipid Class, parsed))
def cls_in(*cs):
    return lambda c, p: c in cs


SETS = [
    # --- membrane structure
    ("Bulk membrane phospholipids (PC, PE)", "Membrane structure",
     "Major structural lipids of cellular membranes", cls_in("PC", "PE")),
    ("Anionic phospholipids (PS, PI, PG, PA, CL, LBPA)", "Membrane structure",
     "Negatively charged lipids; membrane surface charge, protein recruitment", cls_in("PS", "PI", "PG", "PA", "CL", "LBPA")),
    ("Phosphatidylserine (PS)", "Membrane structure",
     "Inner-leaflet lipid; apoptotic 'eat-me' signal, signalling platform", cls_in("PS")),
    ("Phosphatidylinositols (PI, LPI)", "Signalling",
     "Precursors of phosphoinositide signalling (PIPs not detected)", cls_in("PI", "LPI")),
    # --- organelle
    ("Cardiolipin (CL)", "Organelle",
     "Mitochondrial inner membrane; respiration, mitophagy, apoptosis", cls_in("CL")),
    ("Mitochondrial electron transport (CL, CoQ)", "Organelle",
     "Cardiolipin plus coenzyme Q", cls_in("CL", "Co")),
    ("BMP / LBPA (late endosome-lysosome)", "Organelle",
     "Lysobisphosphatidic acid; endolysosomal membranes, lipid sorting", cls_in("LBPA")),
    ("Acylcarnitines (fatty-acid oxidation)", "Organelle",
     "Carnitine shuttle intermediates for mitochondrial beta-oxidation", cls_in("AcCa")),
    # --- storage / sterols
    ("Neutral lipid storage (TG, ChE, ZyE)", "Storage",
     "Lipid droplet core lipids", cls_in("TG", "ChE", "ZyE")),
    ("Triglycerides (TG)", "Storage", "Energy storage; lipid droplets", cls_in("TG")),
    ("Sterol esters (ChE, ZyE)", "Storage",
     "Esterified sterols stored in droplets; cholesterol buffering", cls_in("ChE", "ZyE")),
    ("Cholesterol & sterol homeostasis (Ch, ChE, ZyE)", "Storage",
     "Free and esterified sterols, incl. cholesterol biosynthesis intermediate", cls_in("Ch", "ChE", "ZyE")),
    ("Diacylglycerol & phosphatidic acids (DG, PA, cPA, LPA)", "Signalling",
     "Glycerolipid intermediates and signalling lipids (PKC, mTOR)", cls_in("DG", "PA", "cPA", "LPA")),
    ("Lysophospholipids (LPC, LPE, LPG, LPI, LPA)", "Signalling",
     "PLA2 products / Lands-cycle remodelling intermediates, signalling", cls_in("LPC", "LPE", "LPG", "LPI", "LPA")),
    # --- sphingolipids
    ("Ceramides & sphingoid bases (Cer, SPH)", "Sphingolipids",
     "Stress / apoptosis / inflammatory signalling; sphingolipid hub", cls_in("Cer", "SPH")),
    ("1-deoxy-ceramides (m-base Cer)", "Sphingolipids",
     "Atypical ceramides lacking the 1-OH; non-canonical, cytotoxic", lambda c, p: c == "Cer" and any(pre == "m" for pre, *_ in p["chains"])),
    ("Canonical ceramides (d-base Cer)", "Sphingolipids",
     "Standard dihydroxy sphingoid-base ceramides", lambda c, p: c == "Cer" and any(pre == "d" for pre, *_ in p["chains"])),
    ("Sphingomyelin (SM)", "Sphingolipids",
     "Plasma membrane / raft lipid; source of ceramide", cls_in("SM")),
    ("Glycosphingolipids (HexCer, GM1, GM2, ST)", "Sphingolipids",
     "Glycosylated sphingolipids; rafts, cell recognition, immune signalling", cls_in("Hex1Cer", "GM1", "GM2", "ST")),
    ("Raft-associated lipids (SM, glycosphingolipids, Ch)", "Sphingolipids",
     "Lipids enriched in ordered membrane microdomains", cls_in("SM", "Hex1Cer", "GM1", "GM2", "ST", "Ch")),
    # --- chemical features
    ("Ether lipids (O-/P-)", "Chemical feature",
     "Plasmalogen/ether lipids; peroxisome-derived, antioxidant, ferroptosis", lambda c, p: p["ether"]),
    ("Oxidised / truncated-chain lipids (COOH, CHO)", "Chemical feature",
     "Oxidised phospholipids (annotation-based; may include misannotation)", lambda c, p: p["oxidised"]),
    ("PUFA-containing lipids (any chain >=4 db)", "Chemical feature",
     "Eicosanoid precursors; peroxidation / ferroptosis-prone", lambda c, p: has_pufa(p)),
    ("Arachidonate (20:4)-containing", "Chemical feature",
     "Substrate pool for eicosanoid synthesis", lambda c, p: p["molecular"] and any(x == (20, 4) for x in acyl(p["chains"]))),
    ("DHA (22:6)-containing", "Chemical feature",
     "Omega-3; membrane fluidity, pro-resolving mediator precursor", lambda c, p: p["molecular"] and any(x == (22, 6) for x in acyl(p["chains"]))),
    ("Saturated / monounsaturated acyl lipids (no chain >=2 db)", "Chemical feature",
     "Low unsaturation: stiffer membranes; de novo lipogenesis products",
     lambda c, p: p["molecular"] and c in ("PC", "PE", "PS", "PG", "PI", "TG", "DG", "LBPA", "CL") and bool(acyl(p["chains"])) and all(db <= 1 for _, db in acyl(p["chains"]))),
]
# class-specific PUFA vs non-PUFA splits (molecular-level species only)
for k in ("TG", "PC", "PE", "CL", "PS"):
    SETS.append((f"{k}: PUFA-containing", "Class x unsaturation", f"{k} species with a chain >=4 db",
                 (lambda k: lambda c, p: c == k and has_pufa(p))(k)))
    SETS.append((f"{k}: no PUFA chain", "Class x unsaturation", f"{k} species without a chain >=4 db",
                 (lambda k: lambda c, p: c == k and p["molecular"] and not has_pufa(p))(k)))


def annotate(meta):
    parsed = {i: parse(l) for i, l in meta["Lipid"].items()}
    member = pd.DataFrame({s[0]: [bool(s[3](meta.at[i, "Class"], parsed[i])) for i in meta.index]
                           for s in SETS}, index=meta.index)
    return member


# ---------------------------------------------------------------- stats
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        prev = min(prev, p[i] * n / rank); q[i] = prev
    return q


def set_stats(lfc, member, min_size):
    bg = lfc.median()
    rows = []
    for name, cat, role, _ in SETS:
        inn = lfc[member[name].reindex(lfc.index).fillna(False)]
        out = lfc.drop(inn.index)
        if len(inn) < min_size:
            continue
        p = stats.mannwhitneyu(inn, out, alternative="two-sided").pvalue
        rows.append(dict(set=name, category=cat, role=role, n=len(inn),
                         median_log2FC=inn.median(), relative_log2FC=inn.median() - out.median(),
                         pct_up_2fold=(inn >= 1).mean() * 100, pct_down_2fold=(inn <= -1).mean() * 100,
                         p=p))
    d = pd.DataFrame(rows)
    d["FDR"] = bh(d.p)
    return d


def hit(d):
    """Highlight criterion: FDR<0.01 and |relative log2FC|>=0.5 (lipids are not independent,
    so FDR alone over-calls large sets with tiny shifts)."""
    return (d.FDR < 0.01) & (d.relative_log2FC.abs() >= 0.5)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xlsx")
    ap.add_argument("--outdir", default="lipid_results/sets")
    ap.add_argument("--sheet", default="RT-Averaged Report")
    ap.add_argument("--norm", default="total", choices=["total", "median", "none"])
    ap.add_argument("--min-size", type=int, default=4)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mat, meta = load(args.xlsx, args.sheet)
    norm, _ = normalise(mat, args.norm)
    pseudo = pseudocount(norm)
    member = annotate(meta)

    # annotation tables (editable reference)
    pd.DataFrame([(n, c, r, int(member[n].sum())) for n, c, r, _ in SETS],
                 columns=["set", "category", "biological_role", "n_lipids"]
                 ).to_csv(f"{args.outdir}/set_definitions.csv", index=False)
    long = member.stack()
    long = long[long].reset_index()[["ID", "level_1"]]
    long.columns = ["lipid_ID", "set"]
    long.insert(1, "lipid", meta.loc[long.lipid_ID, "Lipid"].values)
    long.to_csv(f"{args.outdir}/lipid_to_set_annotations.csv", index=False)
    print(f"{len(SETS)} sets defined; {int((member.sum(axis=1) > 0).sum())}/{len(member)} lipids in >=1 set")

    res, lfcs = {}, {}
    for c in ("IFNB", "IFNG"):
        lfc, st = lfc_table(norm, c, "Untreated", pseudo)
        lfc = lfc[st == "both detected"]
        lfcs[c] = lfc
        res[c] = set_stats(lfc, member, args.min_size)
        res[c].sort_values("p").to_csv(f"{args.outdir}/set_results_{c}_vs_Untreated.csv", index=False)
        print(f"\n== {c} vs Untreated: top sets by p ==")
        print(res[c].sort_values("p").head(10)[["set", "n", "relative_log2FC", "p", "FDR"]]
              .round(3).to_string(index=False))

    # ---- dot plot (GO-style)
    order = res["IFNG"].assign(a=lambda d: d.relative_log2FC).sort_values(["category", "a"])
    cats = list(dict.fromkeys(order.category))
    order = pd.concat([order[order.category == k] for k in cats])
    names = order.set.tolist()
    ypos = {n: i for i, n in enumerate(names[::-1])}
    lim = max(res[c].relative_log2FC.abs().max() for c in res) * 1.1
    fig, axs = plt.subplots(1, 2, figsize=(13, 0.33 * len(names) + 2.2), sharey=True)
    for ax, (c, label) in zip(axs, (("IFNB", "IFN-β"), ("IFNG", "IFN-γ"))):
        d = res[c].set_index("set").reindex(names)
        sc = ax.scatter(d.relative_log2FC, [ypos[n] for n in names], s=25 + d.n * 3,
                        c=d.relative_log2FC, cmap="RdBu_r", vmin=-lim, vmax=lim,
                        edgecolor=np.where(hit(d), "k", "#999999"),
                        linewidth=np.where(hit(d), 1.8, .5))
        ax.axvline(0, c="k", lw=.6)
        ax.set_xlim(-lim, lim)
        ax.set_yticks(list(ypos.values())); ax.set_yticklabels(list(ypos.keys()), fontsize=7.5)
        ax.set_xlabel("relative log2FC (set median − all-lipid median)")
        ax.set_title(f"{label} vs Untreated", fontsize=10)
        ax.grid(axis="x", lw=.3)
        ax.set_ylim(-.7, len(names) - .3)
    # category separators
    start = 0
    for k in cats[::-1]:
        n = (order.category == k).sum()
        for ax in axs:
            ax.axhline(start - .5, c="#888", lw=.5)
            ax.text(lim * 1.03, start + n / 2 - .5, k, rotation=270, va="center", fontsize=7, color="#555") if ax is axs[1] else None
        start += n
    fig.colorbar(sc, ax=axs, shrink=.4, label="relative log2FC", pad=.08)
    for n_ in (5, 20, 50):
        axs[1].scatter([], [], s=25 + n_ * 3, c="w", edgecolor="k", label=f"n={n_}")
    axs[1].legend(loc="lower right", fontsize=7, frameon=True, title="set size", title_fontsize=7)
    fig.suptitle("Lipid set analysis (black outline: FDR<0.01 and |rel. log2FC|>=0.5; exploratory, n=1 per condition)", fontsize=10, y=0.93)
    fig.savefig(f"{args.outdir}/lipid_set_dotplot.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- per-set distributions
    fig, axs = plt.subplots(1, 2, figsize=(13, 0.33 * len(names) + 2), sharey=True)
    rng = np.random.default_rng(0)
    for ax, (c, label) in zip(axs, (("IFNB", "IFN-β"), ("IFNG", "IFN-γ"))):
        lfc = lfcs[c]
        for n in names:
            v = lfc[member[n].reindex(lfc.index).fillna(False)]
            y = ypos[n]
            ax.scatter(v, y + rng.uniform(-.25, .25, len(v)), s=5, c="#555", alpha=.5, edgecolor="none")
            ax.plot([v.median()] * 2, [y - .35, y + .35], c="#c0392b", lw=2)
        ax.axvline(0, c="k", lw=.6); ax.axvline(lfc.median(), c="#2471a3", lw=.8, ls="--")
        ax.set_yticks(list(ypos.values())); ax.set_yticklabels(list(ypos.keys()), fontsize=7.5)
        ax.set_xlabel("log2FC vs Untreated"); ax.set_title(f"{label} vs Untreated", fontsize=10)
        ax.set_ylim(-.7, len(names) - .3)
    fig.suptitle("Per-lipid log2FC within each set (red = set median, blue dashed = all-lipid median)", fontsize=10)
    fig.savefig(f"{args.outdir}/lipid_set_distributions.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nWritten to {args.outdir}/")


if __name__ == "__main__":
    main()
