"""Rebuild video/scene.html: inject the processed data and the real figure SVG.

The animation drives the figure's own SVG (output/figure.svg), so every colour,
font and annotation is the figure itself rather than a redrawing of it.
"""

import csv
import re
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VID = ROOT / "video"
S = ["LA1", "LB1", "LC1"]

kos = list(
    csv.DictReader(
        open(ROOT / "data/processed/ko_selected.tsv", encoding="utf-8"), delimiter="\t"
    )
)
classes = list(
    csv.DictReader(
        open(ROOT / "data/processed/class_summary.tsv", encoding="utf-8"),
        delimiter="\t",
    )
)
vals = sorted(float(k[f"rel_abund_{s}"]) for k in kos for s in S)
n = len(vals)

data = {
    "samples": S,
    "vmin": round(vals[int(n * 0.05)], 3),
    "vmax": round(vals[int(n * 0.95) - 1], 3),
    "classes": [
        {
            "code": c["class_code"],
            "name": c["class_name"],
            "group": c["group"],
            "sums": [round(float(c[f"sum_rel_abund_{s}"]), 3) for s in S],
        }
        for c in classes
    ],
    "kos": [
        {
            "id": k["KO"],
            "v": [round(float(k[f"rel_abund_{s}"]), 3) for s in S],
            "d": int(k["n_pathways"]),
            "c": k["classes"].split(";"),
        }
        for k in kos
    ],
}

tpl = (VID / "scene.template.html").read_text(encoding="utf-8")
svg = (ROOT / "output/figure.svg").read_text(encoding="utf-8")
svg = svg[svg.index("<svg") :]
# 根节点只保留 viewBox：pt 单位与 CSS 的 px 混用会把画面裁掉一截
_tag_end = svg.index(">") + 1
_tag = svg[:_tag_end]
_vb = _tag[_tag.index("viewBox="):].split('"')[1]
svg = f'<svg viewBox="{_vb}" preserveAspectRatio="xMidYMin meet" xmlns="http://www.w3.org/2000/svg" version="1.1">' + svg[_tag_end:]

html = tpl.replace("/*__DATA__*/", json.dumps(data, ensure_ascii=False))
html = html.replace("<!--__SVG__-->", svg)
html = html.replace("//<<APPEND>>", "")

out = VID / "scene.html"
out.write_text(html, encoding="utf-8")
print(
    f"{out.name}: {len(html) / 1024:.0f} KB | svg inlined={'<svg' in html} | "
    f"cells={html.count(chr(39) + 'cell-') + html.count(chr(34) + 'cell-')}"
)
