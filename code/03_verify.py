"""Step 3: prove the figure matches the source data. Exits non-zero on any mismatch."""

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
OUT = ROOT / "output"

failures, checks = [], 0


def check(cond, msg):
    global checks
    checks += 1
    if not cond:
        failures.append(msg)
    print(f"{'PASS' if cond else 'FAIL'}  {msg}")


def rows(name):
    return list(csv.DictReader(open(PROC / name, encoding="utf-8"), delimiter="\t"))


kos = rows("ko_selected.tsv")
links = rows("ko_links.tsv")
classes = rows("class_summary.tsv")
prov = json.loads((ROOT / "data" / "provenance.json").read_text(encoding="utf-8"))

# 1. raw inputs are the ones the provenance file claims
for src in prov["sources"]:
    p = ROOT / src["file"]
    check(p.exists(), f"raw input present: {src['file']}")
    if p.exists():
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        check(got == src["sha256"], f"sha256 unchanged: {p.name}")

# 2. every plotted KO exists in both sources
kegg_text = (RAW / "ko00001.keg").read_text(encoding="utf-8", errors="ignore")
src_kos = {
    r[0]
    for r in csv.reader(
        open(RAW / "kegg_KO_profile.xls", encoding="utf-8", errors="ignore"),
        delimiter="\t",
    )
}
check(
    all(re.search(rf"\b{k['KO']}\b", kegg_text) for k in kos),
    "every selected KO appears in KEGG BRITE ko00001",
)
check(
    all(k["KO"] in src_kos for k in kos),
    "every selected KO appears in the figshare abundance table",
)

# 3. tables are internally consistent
ko_ids = {k["KO"] for k in kos}
check(len(ko_ids) == len(kos), "no duplicate KO rows")
check({l["KO"] for l in links} == ko_ids, "link table covers exactly the selected KOs")
for c in classes:
    n = sum(1 for l in links if l["class_code"] == c["class_code"])
    check(int(c["n_links"]) == n, f"link count matches for {c['class_name']}")

for k in kos:
    listed = k["classes"].split(";")
    linked = sorted(l["class_code"] for l in links if l["KO"] == k["KO"])
    check(sorted(listed) == linked, f"classes column matches links for {k['KO']}")

# 4. class sums equal the sum of their member KOs
for c in classes:
    members = [k for k in kos if c["class_code"] in k["classes"].split(";")]
    for s in ("LA1", "LB1", "LC1"):
        want = sum(float(k[f"rel_abund_{s}"]) for k in members)
        got = float(c[f"sum_rel_abund_{s}"])
        check(abs(want - got) < 1e-4, f"sum {c['class_name']} / {s} = {got:.4f}")

# 5. abundances are sane percentages
vals = [float(k[f"rel_abund_{s}"]) for k in kos for s in ("LA1", "LB1", "LC1")]
check(all(0 < v < 100 for v in vals), "all relative abundances in (0, 100)%")
check(len(vals) == 3 * len(kos), "one value per KO per sample")

# 6. figure files were written and are non-trivial
for ext, floor in (("png", 200_000), ("svg", 50_000), ("pdf", 20_000)):
    p = OUT / f"figure.{ext}"
    check(
        p.exists() and p.stat().st_size > floor,
        f"figure.{ext} written ({p.stat().st_size if p.exists() else 0} bytes)",
    )

print(f"\n{checks - len(failures)}/{checks} checks passed")
if failures:
    print("\nFAILURES:")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("all checks passed")
