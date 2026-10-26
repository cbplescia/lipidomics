#!/usr/bin/env python3
"""Lipid process indices, chain landscapes, fatty-acid routing and reaction map.

Contrasts (Mock = Untreated): IFNB vs Mock, IFNG vs Mock, IFNG vs IFNB.

1. Indices / ratios  - enzyme- and process-proxy ratios (desaturation, elongation, PLA2/
                       Lands-cycle remodelling, sphingolipid flux, sterol esterification,
                       neutral-lipid storage, anionic/mitochondrial lipids, oxidation).
2. Chain landscapes  - per class, change in the share of each (total carbons, total double
                       bonds) species within that class (composition, so not affected by
                       how much total lipid a sample has).
3. FA pool / routing - change in each fatty-acyl chain's share within each class, and where
                       20:4, 22:6, 18:1, 16:0 are stored (membrane GP vs neutral lipid ...).
4. Reaction map      - BioPAN-like but descriptive: class nodes coloured by total change,
                       enzymatic edges coloured by change in product/substrate ratio.

Ratios and shares are computed on raw areas within a sample, so overall loading cancels.
Class-to-class ratios are not absolute (ionisation differs) but CHANGES between samples are
comparable. n=1 per condition: descriptive, no p-values.

Usage:
    python lipid_indices.py lipidomics-data.xlsx [--outdir lipid_results/indices]
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch

from lipid_diff_analysis import load
from lipid_set_analysis import parse

SAMPLES = ["Untreated", "IFNB", "IFNG"]
SHOW = {"Untreated": "Mock", "IFNB": "IFN-β", "IFNG": "IFN-γ"}
CONTRASTS = [("IFNB", "Untreated"), ("IFNG", "Untreated"), ("IFNG", "IFNB")]
EPS = 1e-9


def cname(a, b):
    return f"{SHOW[a]} vs {SHOW[b]}"


def log2r(x_a, x_b):
    """log2((a)/(b)) with NaN where either side is zero/NaN"""
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.log2(x_a / x_b)
    return r.replace([np.inf, -np.inf], np.nan) if hasattr(r, "replace") else (r if np.isfinite(r) else np.nan)


# ------------------------------------------------------------------ feature table
def build(path, sheet):
    mat, meta = load(path, sheet)
    mat = mat.fillna(0.0)
    mat.columns = SAMPLES
    rows = []
    for i, r in meta.iterrows():
        p = parse(r["Lipid"])
        ch = p["chains"]
        # fatty-acyl chains usable for pools: unmodified, non-base chains of molecular species,
        # plus the single acyl of mono-acyl classes
        mono = r["Class"] in ("ChE", "LPC", "LPE", "LPI", "LPG", "LPA", "FA")
        acy = [(c, db) for pre, c, db, mod in ch if pre == "" and not mod] if (p["molecular"] or mono) else []
        rows.append(dict(ID=i, Lipid=r["Lipid"], Class=r["Class"], Adduct=r["Adduct"],
                         C=sum(c[1] for c in ch), DB=sum(c[2] for c in ch), has_chain=bool(ch),
                         molecular=p["molecular"], ether=p["ether"], oxidised=p["oxidised"],
                         mbase=any(pre == "m" for pre, *_ in ch), acyl=acy))
    f = pd.DataFrame(rows).set_index("ID")
    return f.join(mat)


def S(f, mask=None, **kw):
    """sum of areas (Series over samples) for rows matching class list / mask"""
    m = pd.Series(True, index=f.index) if mask is None else mask
    if "classes" in kw:
        m &= f.Class.isin(kw["classes"])
    return f.loc[m, SAMPLES].sum()


def pool(f, classes):
    """{chain 'C:db': Series over samples}"""
    out = {}
    for i, r in f[f.Class.isin(classes)].iterrows():
        for c, db in r.acyl:
            k = f"{c}:{db}"
            out[k] = out.get(k, 0) + r[SAMPLES]
    return pd.DataFrame(out).T.astype(float) if out else pd.DataFrame(columns=SAMPLES, dtype=float)


GP = ["PC", "PE", "PS", "PI", "PG", "PA", "LBPA", "CL"]
LYSO = ["LPC", "LPE", "LPI", "LPG", "LPA"]


# ------------------------------------------------------------------ 1. indices
def indices(f):
    P = lambda cl: pool(f, cl)
    gp, tg = P(GP), P(["TG"])
    g = lambda d, k: d.loc[k] if k in d.index else pd.Series(np.nan, index=SAMPLES)
    sat = lambda d: sum(d.loc[k] for k in d.index if k.endswith(":0")) if len(d) else np.nan
    pufa = lambda d: sum(d.loc[k] for k in d.index if int(k.split(":")[1]) >= 4)
    c = lambda *cl: S(f, classes=list(cl))
    mol = f.molecular
    oxm = f.oxidised
    items = [
        # desaturation / elongation
        ("SCD proxy: GP 18:1/18:0", "Desaturation & elongation", g(gp, "18:1") / g(gp, "18:0")),
        ("SCD proxy: GP 16:1/16:0", "Desaturation & elongation", g(gp, "16:1") / g(gp, "16:0")),
        ("SCD proxy: TG 18:1/18:0", "Desaturation & elongation", g(tg, "18:1") / g(tg, "18:0")),
        ("ELOVL6 proxy: GP 18:0/16:0", "Desaturation & elongation", g(gp, "18:0") / g(gp, "16:0")),
        ("FADS/ELOVL: GP 20:4/20:3", "Desaturation & elongation", g(gp, "20:4") / g(gp, "20:3")),
        ("GP 22:6/22:5", "Desaturation & elongation", g(gp, "22:6") / g(gp, "22:5")),
        ("PUFA(>=4db)/saturated, GP chains", "Desaturation & elongation", pufa(gp) / sat(gp)),
        ("PUFA fraction of TG chains", "Desaturation & elongation", pufa(tg) / tg.sum()),
        # remodelling
        ("LPC/PC (PLA2 vs LPCAT)", "Phospholipid remodelling", c("LPC") / c("PC")),
        ("LPE/PE", "Phospholipid remodelling", c("LPE") / c("PE")),
        ("LPI/PI", "Phospholipid remodelling", c("LPI") / c("PI")),
        ("PC/PE", "Phospholipid remodelling", c("PC") / c("PE")),
        ("PS/PE", "Phospholipid remodelling", c("PS") / c("PE")),
        ("PI/PC", "Phospholipid remodelling", c("PI") / c("PC")),
        ("Ether PE / PE (peroxisomal)", "Phospholipid remodelling", S(f, (f.Class == "PE") & f.ether) / c("PE")),
        ("Ether PC / PC", "Phospholipid remodelling", S(f, (f.Class == "PC") & f.ether) / c("PC")),
        ("Oxidised PC / PC", "Phospholipid remodelling", S(f, (f.Class == "PC") & oxm) / c("PC")),
        # sphingolipids
        ("Cer/SM (SMase vs SGMS)", "Sphingolipids", c("Cer") / c("SM")),
        ("SM/PC", "Sphingolipids", c("SM") / c("PC")),
        ("HexCer/Cer (UGCG vs GBA)", "Sphingolipids", c("Hex1Cer") / c("Cer")),
        ("Sphingoid base/Cer (ASAH1)", "Sphingolipids", c("SPH") / c("Cer")),
        ("1-deoxy (m-base) Cer / Cer", "Sphingolipids", S(f, (f.Class == "Cer") & f.mbase) / c("Cer")),
        # sterols
        ("ChE/free Ch (ACAT vs LAL)", "Sterols", c("ChE") / c("Ch")),
        ("Free Ch / (PC+PE+SM)", "Sterols", c("Ch") / c("PC", "PE", "SM")),
        ("ChE / (PC+PE)", "Sterols", c("ChE") / c("PC", "PE")),
        # neutral lipids
        ("TG / (PC+PE)", "Neutral lipid storage", c("TG") / c("PC", "PE")),
        ("DG / TG", "Neutral lipid storage", c("DG") / c("TG")),
        ("DG / PA", "Neutral lipid storage", c("DG") / c("PA")),
        ("PUFA-TG / TG (molecular species)", "Neutral lipid storage",
         S(f, (f.Class == "TG") & f.molecular & f.acyl.map(lambda a: any(db >= 4 for _, db in a))) / S(f, (f.Class == "TG") & f.molecular)),
        ("AcCa / TG", "Neutral lipid storage", c("AcCa") / c("TG")),
        # anionic / mito
        ("CL / PC", "Anionic & mitochondrial", c("CL") / c("PC")),
        ("PG / PC", "Anionic & mitochondrial", c("PG") / c("PC")),
        ("LBPA(BMP), all adducts / PC", "Anionic & mitochondrial", c("LBPA") / c("PC")),
        ("LBPA(BMP) [M-H] only / PC", "Anionic & mitochondrial", S(f, (f.Class == "LBPA") & (f.Adduct == "M-H")) / c("PC")),
        ("PA / PC", "Anionic & mitochondrial", c("PA") / c("PC")),
        ("Anionic GP (PS+PI+PG+PA+CL+LBPA) / total GP", "Anionic & mitochondrial", c("PS", "PI", "PG", "PA", "CL", "LBPA") / c(*GP)),
    ]
    rows = []
    for name, grp, v in items:
        d = dict(index=name, process=grp, **{f"value_{SHOW[s]}": v[s] for s in SAMPLES})
        for a, b in CONTRASTS:
            d[f"log2_{cname(a, b)}"] = np.log2(v[a] / v[b]) if v[a] > 0 and v[b] > 0 else np.nan
        rows.append(d)
    return pd.DataFrame(rows)


def plot_indices(d, outdir):
    cols = [f"log2_{cname(a, b)}" for a, b in CONTRASTS]
    m = d[cols].to_numpy(float)
    lim = np.nanpercentile(np.abs(m), 95) if np.isfinite(m).any() else 1
    lim = max(lim, 0.5)
    fig, ax = plt.subplots(figsize=(7.5, 0.3 * len(d) + 2))
    im = ax.imshow(np.ma.masked_invalid(m), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels([cname(a, b) for a, b in CONTRASTS], fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(d))); ax.set_yticklabels(d["index"], fontsize=7.5)
    for i in range(m.shape[0]):
        for j in range(3):
            if np.isfinite(m[i, j]):
                ax.text(j, i, f"{m[i, j]:+.1f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(m[i, j]) > lim * .65 else "black")
    prev = None
    for i, g in enumerate(d.process):
        if g != prev and i:
            ax.axhline(i - .5, c="k", lw=1)
        prev = g
    gs = list(dict.fromkeys(d.process))
    for g in gs:
        idx = np.where(d.process == g)[0]
        ax.text(3.0, idx.mean(), g, rotation=270, va="center", ha="left", fontsize=7, color="#555")
    fig.colorbar(im, ax=ax, shrink=.4, pad=.22, label="log2 ratio change")
    fig.text(0.5, 0.005, "Descriptive; n=1 per condition. Ratios from raw areas.", ha="center", fontsize=7)
    fig.savefig(f"{outdir}/indices_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ 2. landscapes
def landscape(f, outdir, classes, fname):
    fig, axs = plt.subplots(len(classes), 3, figsize=(11.5, 2.9 * len(classes)), squeeze=False)
    for r, cl in enumerate(classes):
        sub = f[(f.Class == cl) & f.has_chain]
        tot = sub[SAMPLES].sum()
        grid = {s: sub.groupby(["DB", "C"])[s].sum() / tot[s] for s in SAMPLES}
        allk = sorted(set().union(*[g.index for g in grid.values()]))
        if not allk:
            continue
        dbs = sorted({k[0] for k in allk}); cs = sorted({k[1] for k in allk})
        for k, (a, b) in enumerate(CONTRASTS):
            ax = axs[r, k]
            M = np.full((len(dbs), len(cs)), np.nan)
            for (db, c), _ in zip(allk, allk):
                sa, sb = grid[a].get((db, c), 0), grid[b].get((db, c), 0)
                if max(sa, sb) >= 0.005:  # >=0.5% of class in either sample
                    M[dbs.index(db), cs.index(c)] = np.log2((sa + 1e-4) / (sb + 1e-4))
            im = axs[r, k].imshow(np.ma.masked_invalid(M), origin="lower", cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
            ax.set_xticks(range(len(cs))); ax.set_xticklabels(cs, fontsize=6, rotation=90)
            ax.set_yticks(range(len(dbs))); ax.set_yticklabels(dbs, fontsize=6)
            ax.set_title(f"{cl}: {cname(a, b)}", fontsize=8)
            if r == len(classes) - 1:
                ax.set_xlabel("total carbons", fontsize=7)
            if k == 0:
                ax.set_ylabel("total double bonds", fontsize=7)
    fig.colorbar(im, ax=axs, shrink=.3, label="log2 change in share of class")
    fig.suptitle("Chain landscape: change in each (carbons, double bonds) species' share within its class", fontsize=9, y=0.995)
    fig.savefig(f"{outdir}/{fname}", dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ 3. FA pools
def fa_pools(f, outdir):
    cls = ["PC", "PE", "PS", "PI", "PG", "LBPA", "CL", "TG", "DG", "Cer", "ChE"]
    pools = {c: pool(f, [c]) for c in cls}
    allp = pd.concat([p for p in pools.values() if len(p)]).groupby(level=0).sum()
    top = list(allp.sum(axis=1).sort_values(ascending=False).index[:16])
    top = sorted(top, key=lambda k: (int(k.split(":")[0]), int(k.split(":")[1])))
    share = {c: (p.div(p.sum(), axis=1) if len(p) else p) for c, p in pools.items()}
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 0.34 * len(top) + 2.5))
    rows = []
    for ax, (a, b) in zip(axs, CONTRASTS):
        M = np.full((len(top), len(cls)), np.nan)
        for j, c in enumerate(cls):
            s = share[c]
            for i, ch in enumerate(top):
                if ch in s.index:
                    sa, sb = s.loc[ch, a], s.loc[ch, b]
                    if max(sa, sb) >= 0.01:
                        M[i, j] = np.log2((sa + 1e-3) / (sb + 1e-3))
                        rows.append(dict(contrast=cname(a, b), chain=ch, cls=c, share_a=sa, share_b=sb, log2=M[i, j]))
        im = ax.imshow(np.ma.masked_invalid(M), cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
        ax.set_xticks(range(len(cls))); ax.set_xticklabels(cls, fontsize=7, rotation=90)
        ax.set_yticks(range(len(top))); ax.set_yticklabels(top, fontsize=7)
        ax.set_title(cname(a, b), fontsize=9)
    fig.colorbar(im, ax=axs, shrink=.4, label="log2 change in chain's share of class pool")
    fig.suptitle("Fatty-acyl chain composition by lipid class (share of each class's chain pool)", fontsize=10, y=1.0)
    fig.savefig(f"{outdir}/fa_pool_class_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(f"{outdir}/fa_pool_class_changes.csv", index=False)

    # routing of selected chains between compartments (storage vs membrane vs ...)
    pools.update({c: pool(f, [c]) for c in GP + LYSO + ["TG", "DG", "ChE", "Cer"] if c not in pools})
    groups = {"Membrane GP": GP, "Neutral (TG, DG, ChE)": ["TG", "DG", "ChE"], "Lyso": LYSO, "Sphingolipid N-acyl": ["Cer"]}
    gcol = ["#2471a3", "#c0392b", "#16a085", "#8e44ad"]
    sel = ["16:0", "18:1", "20:4", "22:6"]
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.8), sharey=True)
    out = []
    for ax, ch in zip(axs, sel):
        vals = {}
        for g, cl in groups.items():
            tot = pd.Series(0.0, index=SAMPLES)
            for c in cl:
                pc = pools.get(c)
                if pc is not None and ch in pc.index:
                    tot = tot + pc.loc[ch]
            vals[g] = tot
        df = pd.DataFrame(vals).T
        frac = df.div(df.sum()) * 100
        bottom = np.zeros(3)
        for (g, row), col in zip(frac.iterrows(), gcol):
            ax.bar([SHOW[s] for s in SAMPLES], row.values, bottom=bottom, color=col, label=g)
            bottom += row.values
        ax.set_title(f"{ch} chain", fontsize=9)
        for g, row in frac.iterrows():
            out.append(dict(chain=ch, group=g, **{SHOW[s]: row[s] for s in SAMPLES}))
    axs[0].set_ylabel("% of that chain's total signal"); axs[-1].legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.suptitle("Where each fatty-acyl chain is stored (class groups; ionisation differs by class, compare across samples)", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{outdir}/fa_routing.png", dpi=200); plt.close(fig)
    pd.DataFrame(out).to_csv(f"{outdir}/fa_routing.csv", index=False)


# ------------------------------------------------------------------ 4. reaction map
NODES = {
    "PC": (0, 0), "LPC": (-2.6, 1.6), "PA": (2.6, 1.6), "DG": (5.2, 0), "TG": (7.8, 0),
    "PE": (0, -3), "PS": (-2.8, -4.4), "SM": (-2.6, 3.8), "Cer": (0, 5.4), "Hex1Cer": (2.8, 5.4),
    "SPH": (-2.6, 6.8), "Ch": (-6, -1.2), "ChE": (-8.6, -1.2), "PG": (3.2, -3), "CL": (6, -3),
}
EDGES = [  # substrate, product, enzymes
    ("PC", "LPC", "PLA2"), ("LPC", "PC", "LPCAT"), ("PC", "PA", "PLD"), ("PA", "DG", "LPIN/PLPP"),
    ("DG", "PA", "DGK"), ("DG", "TG", "DGAT"), ("TG", "DG", "PNPLA2/LIPA"), ("DG", "PC", "CHPT/CEPT"),
    ("PC", "SM", "SGMS"), ("SM", "Cer", "SMPD"), ("Cer", "Hex1Cer", "UGCG"), ("Hex1Cer", "Cer", "GBA1"),
    ("Cer", "SPH", "ASAH1"), ("SPH", "Cer", "CERS"), ("PS", "PE", "PISD"), ("PE", "PC", "PEMT"),
    ("PC", "PS", "PTDSS1"), ("Ch", "ChE", "SOAT1"), ("ChE", "Ch", "LIPA/NCEH1"), ("PG", "CL", "CRLS1"),
]


def reaction_map(f, outdir):
    SC = 1.55
    pos = {n: (x * SC, y * SC) for n, (x, y) in NODES.items()}
    cs = {c: S(f, classes=[c]) for c in NODES}
    tot = f[SAMPLES].sum()
    fig, axs = plt.subplots(3, 1, figsize=(15, 29))
    cmap = plt.get_cmap("RdBu_r")
    pairs = {(e[0], e[1]) for e in EDGES}
    rows = []
    for ax, (a, b) in zip(axs, CONTRASTS):
        for s, p, enz in EDGES:
            if min(cs[s][a], cs[s][b], cs[p][a], cs[p][b]) <= 0:
                continue
            v = np.log2((cs[p][a] / cs[s][a]) / (cs[p][b] / cs[s][b]))
            rows.append(dict(contrast=cname(a, b), substrate=s, product=p, enzymes=enz, log2_ratio_change=v))
            col = cmap((np.clip(v, -3, 3) + 3) / 6)
            (x1, y1), (x2, y2) = pos[s], pos[p]
            rev = (p, s) in pairs
            rad = 0.28 if rev else 0.0
            ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=16, lw=3, color=col,
                                         connectionstyle=f"arc3,rad={rad}", shrinkA=26, shrinkB=26, zorder=2))
            dx, dy = x2 - x1, y2 - y1
            mx, my = (x1 + x2) / 2 + (0.5 * rad + 0.06 * (1 if rev else 0)) * dy, (y1 + y2) / 2 - (0.5 * rad + 0.06 * (1 if rev else 0)) * dx
            if not rev:
                L = max((dx ** 2 + dy ** 2) ** .5, 1e-6)
                mx, my = mx + dy / L * 0.55, my - dx / L * 0.55
            ax.text(mx, my, f"{enz}\n{v:+.1f}", fontsize=7.5, ha="center", va="center", color="#222", zorder=5,
                    bbox=dict(fc="white", ec=col, lw=1, alpha=.9, boxstyle="round,pad=0.15"))
        for n, (x, y) in pos.items():
            v = np.log2((cs[n][a] / tot[a]) / (cs[n][b] / tot[b])) if cs[n][a] > 0 and cs[n][b] > 0 else np.nan
            col = cmap((np.clip(v, -3, 3) + 3) / 6) if np.isfinite(v) else "#dddddd"
            ax.scatter([x], [y], s=1500, c=[col], edgecolor="k", zorder=3)
            ax.text(x, y, f"{n}\n{v:+.1f}" if np.isfinite(v) else n, ha="center", va="center", fontsize=8.5, zorder=4, fontweight="bold")
        ax.set_title(cname(a, b), fontsize=13)
        ax.set_xlim(-15, 14.5); ax.set_ylim(-9, 11.8); ax.set_aspect("equal"); ax.axis("off")
    fig.suptitle("Descriptive reaction map (n=1)\nnode = class share of total signal, log2 change; "
                 "edge = change in product/substrate ratio (red = shifted toward product, blue = toward substrate)",
                 fontsize=11, y=0.905)
    fig.savefig(f"{outdir}/reaction_map.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(f"{outdir}/reaction_edges.csv", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xlsx")
    ap.add_argument("--outdir", default="lipid_results/indices")
    ap.add_argument("--sheet", default="RT-Averaged Report")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    f = build(args.xlsx, args.sheet)
    d = indices(f)
    d.round(4).to_csv(f"{args.outdir}/indices.csv", index=False)
    plot_indices(d, args.outdir)
    landscape(f, args.outdir, ["PC", "PE", "PS", "LBPA"], "chain_landscape_GP.png")
    landscape(f, args.outdir, ["TG", "CL", "SM", "Cer"], "chain_landscape_TG_CL_SP.png")
    fa_pools(f, args.outdir)
    reaction_map(f, args.outdir)
    print(f"Written to {args.outdir}/")


if __name__ == "__main__":
    main()
