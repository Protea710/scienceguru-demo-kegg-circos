"""Rebuild video/scene.html from the processed tables and the HTML template.

The animation must show the same numbers as the figure, so the data is injected
from data/processed/ rather than typed by hand.
"""

import csv
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
html = tpl.replace("/*__DATA__*/", json.dumps(data, ensure_ascii=False)).replace(
    "//<<APPEND>>", ""
)
out = VID / "scene.html"
out.write_text(html, encoding="utf-8")
print(
    f"{out.name}: {len(html) / 1024:.0f} KB, {len(data['kos'])} KOs, {len(data['classes'])} classes"
)
