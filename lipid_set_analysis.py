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

# ---- process / compartment / biophysics collection -------------------------------
# Class-level assignments from the lipid-biology literature. Whole-cell lipidomics
# cannot resolve leaflet or organelle, so these are "lipids that characteristically
# live there / do that", not measurements of it.
PL = ("PC", "PE", "PS", "PI", "PG", "PA", "CL", "LBPA", "LPC", "LPE", "LPI", "LPG", "LPA")
PUFA_FA = {(20, 4), (20, 5), (22, 4), (22, 5), (22, 6)}


def has_fa(p, wanted):
    return p["molecular"] and any(x in wanted for x in acyl(p["chains"]))


SETS_PROCESS = [
    # signalling processes
    ("Apoptosis / efferocytosis signalling (PS, Cer, SPH, CL)", "Signalling process",
     "PS exposure ('eat-me'), ceramide/sphingosine pro-apoptotic, CL externalisation", cls_in("PS", "Cer", "SPH", "CL")),
    ("Mitophagy & mitochondrial apoptosis (CL, PG, CoQ)", "Signalling process",
     "Cardiolipin as mitophagy/apoptosis platform; PG = CL precursor", cls_in("CL", "PG", "Co")),
    ("PKC / mTOR second messengers (DG, PA)", "Signalling process",
     "Diacylglycerol and phosphatidic acid; PKC, mTORC1 signalling", cls_in("DG", "PA")),
    ("Lysophospholipid GPCR signalling (LPA, cPA, LPC, LPI)", "Signalling process",
     "Extracellular lysolipid mediators acting on LPAR/G2A/GPR55; chemotaxis, inflammation", cls_in("LPA", "cPA", "LPC", "LPI")),
    ("Sphingolipid stress response (Cer, SPH)", "Signalling process",
     "Ceramide / sphingoid bases in stress, cytokine and apoptosis signalling", cls_in("Cer", "SPH")),
    ("Eicosanoid / pro-resolving mediator precursor pool (20:4, 20:5, 22:4-6 phospholipids)", "Signalling process",
     "PUFA in membrane phospholipids released by PLA2 to make eicosanoids / resolvins",
     lambda c, p: c in PL and has_fa(p, PUFA_FA)),
    ("PUFA stored in neutral lipids (TG, DG, ChE with 20:4-22:6)", "Signalling process",
     "Lipid-droplet PUFA reservoir (ATGL/lipophagy) feeding mediator synthesis",
     lambda c, p: c in ("TG", "DG", "ChE") and any(x in PUFA_FA for x in acyl(p["chains"]))),
    ("Ferroptosis-relevant lipids (AA/AdA-PE, PUFA ether PE/PC)", "Signalling process",
     "PUFA-PE and PUFA ether lipids are the peroxidation substrates in ferroptosis",
     lambda c, p: (c == "PE" and has_fa(p, {(20, 4), (22, 4)})) or (p["ether"] and c in ("PE", "PC") and total_db(p) >= 4)),
    # compartments / leaflets
    ("Cytosolic (inner) leaflet: PM & endosomes (PS, PE, PI, PA)", "Compartment / leaflet",
     "Aminophospholipids and PI/PA are enriched on the cytosolic face of plasma and endosomal membranes", cls_in("PS", "PE", "PI", "PA")),
    ("Exoplasmic / lumenal leaflet: PM outer & endosome lumen (PC, SM, glycosphingolipids)", "Compartment / leaflet",
     "Choline lipids and sphingolipids on the outer / lumenal face", cls_in("PC", "SM", "Hex1Cer", "GM1", "GM2", "ST")),
    ("Late endosome / ILV / lysosome (LBPA, Cer, Ch, GM1, GM2)", "Compartment / leaflet",
     "BMP/LBPA-rich intralumenal vesicles; ceramide and cholesterol sorting; gangliosides", cls_in("LBPA", "Cer", "Ch", "GM1", "GM2")),
    ("Mitochondrial membranes (CL, PG, CoQ)", "Compartment / leaflet",
     "Inner-membrane lipids: cardiolipin, PG (CL precursor), ubiquinone", cls_in("CL", "PG", "Co")),
    ("Ordered raft-like domains (SM, glycosphingolipids, Ch)", "Compartment / leaflet",
     "Saturated sphingolipid/cholesterol microdomains; receptor signalling platforms", cls_in("SM", "Hex1Cer", "GM1", "GM2", "ST", "Ch")),
    # biophysics
    ("Negative curvature / cone-shaped (PE, PA, DG, CL, Cer, Ch)", "Membrane biophysics",
     "Promote fusion, fission, tubulation and non-bilayer intermediates", cls_in("PE", "PA", "DG", "CL", "Cer", "Ch")),
    ("Positive curvature / inverted-cone (lysophospholipids, gangliosides)", "Membrane biophysics",
     "Promote budding, micelle-like and highly curved structures", cls_in("LPC", "LPE", "LPI", "LPG", "LPA", "GM1", "GM2", "ST")),
    ("Bilayer-forming cylindrical lipids (PC, PS, PI, PG, SM)", "Membrane biophysics",
     "Planar bilayer scaffolds", cls_in("PC", "PS", "PI", "PG", "SM")),
    ("Fluidising: polyunsaturated phospholipids (chain >=4 db)", "Membrane biophysics",
     "Increase disorder, flexibility and permeability", lambda c, p: c in PL and has_pufa(p)),
    ("Rigidifying: saturated/MUFA phospholipids (no chain >=2 db)", "Membrane biophysics",
     "Tighter packing, higher melting temperature",
     lambda c, p: c in PL and p["molecular"] and bool(acyl(p["chains"])) and all(db <= 1 for _, db in acyl(p["chains"]))),
    ("Ordering / condensing lipids (Ch, SM, Cer, glycosphingolipids)", "Membrane biophysics",
     "Increase packing and thickness, reduce fluidity", cls_in("Ch", "SM", "Cer", "Hex1Cer", "GM1", "GM2", "ST")),
    ("Long acyl chain (>=C22) lipids: thickness / interdigitation", "Membrane biophysics",
     "Very-long chains span or interdigitate across leaflets",
     lambda c, p: p["molecular"] and bool(acyl(p["chains"])) and any(cc >= 22 for cc, _ in acyl(p["chains"])) and c not in ("TG", "ChE")),
    ("Short / medium acyl chain (<C16) lipids: increased fluidity", "Membrane biophysics",
     "Short chains thin and loosen membranes; also fatty-acid oxidation intermediates",
     lambda c, p: p["molecular"] and bool(acyl(p["chains"])) and any(cc < 16 for cc, _ in acyl(p["chains"])) and c not in ("AcCa",)),
    ("Ether-linked phospholipids: fusion-prone, tighter packing", "Membrane biophysics",
     "Plasmalogens favour hexagonal phases and fusion; reduce bilayer permeability", lambda c, p: p["ether"] and c in PL),
]
SETS += SETS_PROCESS

# ---- subcellular localization collection ---------------------------------------
# Lipids characteristically ENRICHED in (or made in / marking) each compartment.
# Whole-cell lipidomics cannot measure localization; most lipids occur in several
# compartments, so these are marker/enrichment-based hypotheses, strongest for
# cardiolipin (mitochondria), LBPA/BMP (late endosome/lysosome), TG/ChE (droplets).
LOC = {
    "mito": "Mitochondrial inner membrane (CL, PG, CoQ)",
    "fao": "Mitochondrial fatty-acid oxidation (acylcarnitines)",
    "endo": "Endolysosome (LBPA/BMP, GM1, GM2, SPH)",
    "ld": "Lipid droplets (TG, ChE, ZyE core)",
    "pm_out": "Plasma membrane outer leaflet / rafts (SM, glycosphingolipids, Ch)",
    "pm_in": "Plasma membrane inner leaflet (PS, PI, PA)",
    "er": "ER lipid synthesis intermediates (DG, PA, PI, ZyE)",
    "golgi": "Golgi sphingolipid maturation (SM, HexCer, ST, GM1/GM2)",
    "perox": "Peroxisome-dependent ether lipids (O-/P-)",
    "extra": "Extracellular / secreted lysolipid mediators (LPA, cPA, LPC, LPI)",
}
SETS_LOC = [
    (LOC["mito"], "Subcellular localization", "Cardiolipin is the hallmark mitochondrial inner-membrane lipid; PG is its precursor; CoQ shuttles electrons", cls_in("CL", "PG", "Co")),
    (LOC["fao"], "Subcellular localization", "Carnitine shuttle moves fatty acids into mitochondria for beta-oxidation", cls_in("AcCa")),
    (LOC["endo"], "Subcellular localization", "BMP/LBPA marks late endosome ILVs; gangliosides and sphingosine are lysosomal degradation intermediates", cls_in("LBPA", "GM1", "GM2", "SPH")),
    (LOC["ld"], "Subcellular localization", "Neutral-lipid core of lipid droplets", cls_in("TG", "ChE", "ZyE")),
    (LOC["pm_out"], "Subcellular localization", "Outer-leaflet sphingolipids and cholesterol of the plasma membrane", cls_in("SM", "Hex1Cer", "GM1", "GM2", "ST", "Ch")),
    (LOC["pm_in"], "Subcellular localization", "Anionic lipids on the cytosolic face of the plasma membrane", cls_in("PS", "PI", "PA")),
    (LOC["er"], "Subcellular localization", "Glycerolipid / sterol biosynthetic intermediates made in the ER", cls_in("DG", "PA", "PI", "ZyE")),
    (LOC["golgi"], "Subcellular localization", "SM, glucosylceramide and complex GSL are made in the Golgi", cls_in("SM", "Hex1Cer", "ST", "GM1", "GM2")),
    (LOC["perox"], "Subcellular localization", "Ether-lipid synthesis starts in peroxisomes", lambda c, p: p["ether"]),
    (LOC["extra"], "Subcellular localization", "Lysolipid mediators act outside the cell on GPCRs", cls_in("LPA", "cPA", "LPC", "LPI")),
]
SETS += SETS_LOC

COLLECTIONS = {
    "functional": ["Membrane structure", "Organelle", "Storage", "Signalling", "Sphingolipids",
                   "Chemical feature", "Class x unsaturation"],
    "process_biophysics": ["Signalling process", "Compartment / leaflet", "Membrane biophysics"],
    "localization": ["Subcellular localization"],
}


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


# ---------------------------------------------------------------- plots
def make_plots(tag, cats_in, res, lfcs, member, outdir):
    sub = {c: res[c][res[c].category.isin(cats_in)] for c in res}
    order = sub["IFNG"].assign(a=lambda d: d.relative_log2FC).sort_values(["category", "a"])
    cats = [k for k in cats_in if k in set(order.category)]
    order = pd.concat([order[order.category == k] for k in cats])
    names = order.set.tolist()
    ypos = {n: i for i, n in enumerate(names[::-1])}
    lim = max(sub[c].relative_log2FC.abs().max() for c in sub) * 1.1
    fig, axs = plt.subplots(1, 2, figsize=(15, 0.36 * len(names) + 2.2), sharey=True)
    for ax, (c, label) in zip(axs, (("IFNB", "IFN-β"), ("IFNG", "IFN-γ"))):
        d = sub[c].set_index("set").reindex(names)
        sc = ax.scatter(d.relative_log2FC, [ypos[n] for n in names], s=25 + d.n * 3,
                        c=d.relative_log2FC, cmap="RdBu_r", vmin=-lim, vmax=lim,
                        edgecolor=np.where(hit(d), "k", "#999999"),
                        linewidth=np.where(hit(d), 1.8, .5))
        ax.axvline(0, c="k", lw=.6)
        ax.set_xlim(-lim, lim)
        ax.set_yticks(list(ypos.values())); ax.set_yticklabels(list(ypos.keys()), fontsize=7)
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
    fig.savefig(f"{outdir}/{tag}_dotplot.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


    # ---- per-set distributions
    fig, axs = plt.subplots(1, 2, figsize=(15, 0.36 * len(names) + 2), sharey=True)
    rng = np.random.default_rng(0)
    for ax, (c, label) in zip(axs, (("IFNB", "IFN-β"), ("IFNG", "IFN-γ"))):
        lfc = lfcs[c]
        for n in names:
            v = lfc[member[n].reindex(lfc.index).fillna(False)]
            y = ypos[n]
            ax.scatter(v, y + rng.uniform(-.25, .25, len(v)), s=5, c="#555", alpha=.5, edgecolor="none")
            ax.plot([v.median()] * 2, [y - .35, y + .35], c="#c0392b", lw=2)
        ax.axvline(0, c="k", lw=.6); ax.axvline(lfc.median(), c="#2471a3", lw=.8, ls="--")
        ax.set_yticks(list(ypos.values())); ax.set_yticklabels(list(ypos.keys()), fontsize=7)
        ax.set_xlabel("log2FC vs Untreated"); ax.set_title(f"{label} vs Untreated", fontsize=10)
        ax.set_ylim(-.7, len(names) - .3)
    fig.suptitle("Per-lipid log2FC within each set (red = set median, blue dashed = all-lipid median)", fontsize=10)
    fig.savefig(f"{outdir}/{tag}_distributions.png", dpi=200, bbox_inches="tight")
    plt.close(fig)



def make_cell_schematic(res, outdir):
    """Cartoon cell with each compartment coloured by its set's relative log2FC."""
    from matplotlib.patches import Ellipse, Circle, Arc
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    sub = {c: res[c].set_index("set") for c in res}
    lim = 2.5  # colour scale saturates here so modest shifts stay visible
    cmap, norm_ = plt.get_cmap("RdBu_r"), Normalize(-lim, lim)

    fig, axs = plt.subplots(2, 1, figsize=(13, 12.5))
    for ax, (c, label) in zip(axs, (("IFNB", "IFN-β"), ("IFNG", "IFN-γ"))):
        def col(key):
            r = sub[c].relative_log2FC.get(LOC[key], np.nan)
            ok = bool((sub[c].FDR.get(LOC[key], 1) < 0.01) and abs(r) >= 0.5)
            return (cmap(norm_(r)) if np.isfinite(r) else "#dddddd"), r, ok

        def tag(key, anchor, text_xy, ha):
            _, r, ok = col(key)
            ax.annotate(f"{LOC[key].split(' (')[0]}\n{r:+.2f}{' *' if ok else ''}", anchor, text_xy,
                        fontsize=9, ha=ha, va="center", fontweight="bold" if ok else "normal",
                        arrowprops=dict(arrowstyle="-", color="#888", lw=.8, shrinkA=2, shrinkB=0))

        ec = lambda ok: "k" if ok else "#666666"
        f, r, ok = col("pm_out")
        ax.add_patch(Ellipse((0, 0), 11.4, 8.6, fc="#fbfbfb", ec=f, lw=9))
        f, r, ok = col("pm_in")
        ax.add_patch(Ellipse((0, 0), 10.4, 7.6, fc="none", ec=f, lw=5))
        ax.add_patch(Circle((-2.3, 0.6), 1.35, fc="#e7e7ee", ec="#999999", lw=1))
        ax.text(-2.3, 0.6, "Nucleus\n(not assessed)", ha="center", va="center", fontsize=7, color="#777")
        f, r, ok = col("er")
        ax.add_patch(Ellipse((-2.3, 0.6), 4.3, 4.1, fc="none", ec=f, lw=7, ls=(0, (1, 0.5))))
        f, r, ok = col("golgi")
        for k in range(4):
            ax.add_patch(Arc((1.4, 2.0 - 0.28 * k), 2.6 - 0.2 * k, 1.1, theta1=200, theta2=340, ec=f, lw=6, capstyle="round"))
        f, r, ok = col("mito"); ff, rf, okf = col("fao")
        for (x, y, a) in ((2.5, -1.7, 20), (-0.6, -2.7, -10), (3.4, 0.1, 80)):
            ax.add_patch(Ellipse((x, y), 2.0, 0.95, angle=a, fc=f, ec=ec(ok), lw=1.5))
            ax.add_patch(Ellipse((x, y), 1.1, 0.4, angle=a, fc=ff, ec=ec(okf), lw=1))
        f, r, ok = col("endo")
        for (x, y) in ((0.3, -0.3), (1.3, -0.9), (0.2, 1.0)):
            ax.add_patch(Circle((x, y), 0.4, fc=f, ec=ec(ok), lw=1.5))
        f, r, ok = col("ld")
        for (x, y, rad) in ((-4.0, -2.0, .45), (-3.3, -2.4, .33), (-3.6, -1.4, .3), (-4.3, -1.3, .22)):
            ax.add_patch(Circle((x, y), rad, fc=f, ec=ec(ok), lw=1.5))
        f, r, ok = col("perox")
        for (x, y) in ((4.0, 2.0), (-0.3, 3.1), (-1.5, -2.2)):
            ax.add_patch(Circle((x, y), 0.25, fc=f, ec=ec(ok), lw=1.2))
        f, r, ok = col("extra")
        for (x, y) in ((6.5, 1.0), (7.0, 1.6), (6.6, 2.2), (7.3, 2.6)):
            ax.add_patch(Circle((x, y), 0.17, fc=f, ec=ec(ok), lw=1))

        tag("pm_out", (4.0, 3.1), (9.6, 4.6), "left")
        tag("perox", (4.0, 2.0), (9.6, 3.6), "left")
        tag("extra", (7.3, 2.0), (9.6, 2.3), "left")
        tag("mito", (3.8, 0.3), (9.6, 0.0), "left")
        tag("fao", (2.7, -1.9), (9.6, -1.5), "left")
        tag("golgi", (1.4, 1.9), (1.4, 5.7), "center")
        tag("er", (-4.4, 1.5), (-7.8, 2.2), "right")
        tag("ld", (-4.3, -1.9), (-7.8, -1.5), "right")
        tag("pm_in", (-3.8, -2.8), (-7.8, -3.9), "right")
        tag("endo", (1.3, -1.3), (2.2, -5.8), "center")
        ax.set_xlim(-12.5, 13.5); ax.set_ylim(-6.6, 6.6); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"{label} vs Untreated", fontsize=12)
    cax = fig.add_axes([0.97, 0.35, 0.015, 0.3])
    cb = fig.colorbar(ScalarMappable(norm_, cmap), cax=cax, extend="both")
    cb.set_label("relative log2FC (set median − all-lipid median)")
    fig.suptitle("Lipid changes by subcellular localization (class-level assignment; * = FDR<0.01 & |rel. log2FC|≥0.5; exploratory, n=1)",
                 fontsize=10, y=0.92)
    fig.savefig(f"{outdir}/localization_cell_schematic.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


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

    for tag, cats_in in COLLECTIONS.items():
        make_plots(tag, cats_in, res, lfcs, member, args.outdir)
    make_cell_schematic(res, args.outdir)
    print(f"\nWritten to {args.outdir}/")


if __name__ == "__main__":
    main()
