"""Rendering engine for lipid lifecycle pathway diagrams (synthesis -> remodeling -> degradation -> recycling)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

COLS = [
    ("1  SYNTHESIS", "#DCEBFA"),
    ("2  REMODELING / CONVERSION", "#FDF0D5"),
    ("3  DEGRADATION", "#FADBD8"),
    ("4  RECYCLING / FATE", "#DDF2DD"),
]
COLW, ROWH = 5.6, 1.9
STYLE = {
    "self": dict(fc="#1F4E79", ec="#1F4E79", tc="white", weight="bold", ls="-"),
    "met": dict(fc="#F5F5F5", ec="#555555", tc="#222222", weight="normal", ls="-"),
    "rec": dict(fc="#FFFFFF", ec="#2E7D32", tc="#1B5E20", weight="normal", ls="--"),
    "ext": dict(fc="#FFFFFF", ec="#777777", tc="#333333", weight="normal", ls=":"),
}


def render(path, title, subtitle, nodes, edges, note=""):
    """nodes: {id: (col,row,text,kind)}; edges: [(a,b,label,rad,kind)] kind 'n' normal / 'r' recycling."""
    nrows = max(n[1] for n in nodes.values()) + 1
    W = COLW * 4
    H = nrows * ROWH + 2.3
    fig = plt.figure(figsize=(W + 0.6, H), dpi=100)
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
    ax.set_xlim(-0.2, W + 0.2)
    ax.set_ylim(-(nrows * ROWH + 0.9), 1.6)
    ax.set_aspect("equal")
    ax.axis("off")
    bottom = -(nrows * ROWH + 0.5)
    for i, (name, color) in enumerate(COLS):
        ax.add_patch(Rectangle((i * COLW + 0.05, bottom), COLW - 0.1, -bottom + 0.35,
                               fc=color, ec="none", zorder=0))
        ax.text(i * COLW + COLW / 2, 0.15, name, ha="center", va="center",
                fontsize=12, weight="bold", color="#333333")
    ax.text(0.1, 1.35, title, fontsize=22, weight="bold", color="#1F4E79", va="center")
    ax.text(0.1, 0.85, subtitle, fontsize=11.5, color="#444444", va="center", style="italic")

    used = {e[0] for e in edges} | {e[1] for e in edges}
    nodes = {k: v for k, v in nodes.items() if k in used}
    pos, txt = {}, {}
    for k, (c, r, text, kind) in nodes.items():
        s = STYLE[kind]
        x, y = c * COLW + COLW / 2, -(r * ROWH + 1.0)
        pos[k] = (x, y)
        txt[k] = ax.text(x, y, text, ha="center", va="center", fontsize=12.5, weight=s["weight"],
                         color=s["tc"], zorder=3, linespacing=1.25,
                         bbox=dict(boxstyle="round,pad=0.45", fc=s["fc"], ec=s["ec"], lw=1.6, ls=s["ls"]))
    fig.canvas.draw()
    for a, b, label, rad, kind in edges:
        (x1, y1), (x2, y2) = pos[a], pos[b]
        col = "#2E7D32" if kind == "r" else "#444444"
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), zorder=1,
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.7, mutation_scale=16,
                                    ls="--" if kind == "r" else "-", shrinkA=3, shrinkB=3,
                                    patchA=txt[a].get_bbox_patch(), patchB=txt[b].get_bbox_patch(),
                                    connectionstyle=f"arc3,rad={rad}"))
        if label:
            mx = (x1 + x2) / 2 + 0.5 * rad * (y2 - y1)
            my = (y1 + y2) / 2 - 0.5 * rad * (x2 - x1)
            ax.text(mx, my, label, ha="center", va="center", fontsize=10.5, style="italic",
                    color=col, zorder=4, linespacing=1.15,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.92))
    foot = "Solid arrows: metabolic conversion (enzymes in italics).  Dashed green: recycling / salvage.  Simplified schematic."
    if note:
        foot = note + "\n" + foot
    ax.text(0.1, bottom - 0.3, foot, fontsize=9, color="#555555", va="top")
    fig.savefig(path, dpi=100, facecolor="white")
    plt.close(fig)
