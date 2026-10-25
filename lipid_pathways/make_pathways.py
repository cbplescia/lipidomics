"""Generate lifecycle pathway PNGs: python make_pathways.py [out_dir]"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathway_engine import render
from pathways import P, build

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "png")
os.makedirs(out, exist_ok=True)
for name, d in P.items():
    nodes, edges = build(d["nodes"], d["edges"])
    render(os.path.join(out, f"{name}_lifecycle.png"), d["title"], d["sub"], nodes, edges)
    print("wrote", name)
