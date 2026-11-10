"""Step 1: fetch the two public sources and build the three tables used by the plot.

Sources
-------
1. KEGG BRITE ko00001 (KO -> pathway class hierarchy), rest.kegg.jp
2. KEGG ortholog abundance profile, figshare 10.6084/m9.figshare.29645648 (CC BY 4.0)

Everything written here is derived from those two files; nothing is simulated.
"""

import csv
import hashlib
import json
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"

KEGG_URL = "https://rest.kegg.jp/get/br:ko00001"
FIGSHARE_FILE = "https://ndownloader.figshare.com/files/56550593"
FIGSHARE_DOI = "10.6084/m9.figshare.29645648.v1"

# The 17 level-2 categories drawn on the outer arc of the reference figure.
CLASSES = [
    ("09104", "Nucleotide metabolism"),
    ("09102", "Energy metabolism"),
    ("09101", "Carbohydrate metabolism"),
    ("09106", "Metabolism of other amino acids"),
    ("09105", "Amino acid metabolism"),
    ("09111", "Xenobiotics biodegradation"),
    ("09103", "Lipid metabolism"),
    ("09122", "Translation"),
    ("09124", "Replication and repair"),
    ("09123", "Folding, sorting and degradation"),
    ("09131", "Membrane transport"),
    ("09132", "Signal transduction"),
    ("09143", "Cell growth and death"),
    ("09145", "Cellular community - prokaryotes"),
    ("09149", "Aging"),
    ("09171", "Infectious disease: bacterial"),
    ("09175", "Drug resistance: antimicrobial"),
]
KEEP = {c for c, _ in CLASSES}
METABOLISM = {"09104", "09102", "09101", "09106", "09105", "09111", "09103"}

# One soil site (L), three amendment treatments: A control, B manure, C straw.
SAMPLES = ["LA1", "LB1", "LC1"]
MAX_PER_CLASS = 4
MAX_KO = 52

PATH_RE = re.compile(r"^C\s+(\d{5})\s")


def sha256(path):
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def download(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=180) as r, open(dest, "wb") as f:
        f.write(r.read())


def parse_kegg(path):
    """Return ko -> dict(classes=[...], pathways=[...]) and pathway -> class."""
    ko2cl = defaultdict(list)
    ko2path = defaultdict(set)
    path_cls = {}
    cur_cls = cur_path = None
    for line in open(path, encoding="utf-8", errors="ignore"):
        m_path = PATH_RE.match(line)
        if m_path:
            cur_path = m_path.group(1)
            if cur_cls:
                path_cls[cur_path] = cur_cls
        elif line.startswith("B") and line[1] == " ":
            cur_cls = line[3:8].strip()
            cur_path = None
        elif line.startswith("D"):
            m = re.match(r"D\s+(K\d{5})", line)
            if m:
                ko = m.group(1)
                if cur_path:
                    ko2path[ko].add(cur_path)
                if cur_cls and cur_cls not in ko2cl[ko]:
                    ko2cl[ko].append(cur_cls)
    return ko2cl, ko2path, path_cls


def load_matrix(path):
    rows = list(
        csv.reader(open(path, encoding="utf-8", errors="ignore"), delimiter="\t")
    )
    header = rows[0]
    samples = [h for h in header if re.fullmatch(r"[A-Z]{2}\d", h)]
    pos = {h: i for i, h in enumerate(header)}
    mat = {}
    for r in rows[1:]:
        if not r or not r[0].startswith("K"):
            continue
        try:
            mat[r[0]] = {s: float(r[pos[s]] or 0) for s in samples}
        except (ValueError, IndexError):
            continue
    return mat, samples


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    PROC.mkdir(parents=True, exist_ok=True)

    kegg_file = RAW / "ko00001.keg"
    ab_file = RAW / "kegg_KO_profile.xls"
    print(f"fetching KEGG BRITE  -> {kegg_file.name}")
    download(KEGG_URL, kegg_file)
    print(f"fetching figshare    -> {ab_file.name}")
    download(FIGSHARE_FILE, ab_file)

    ko2cl, ko2path, path_cls = parse_kegg(kegg_file)
    mat, all_samples = load_matrix(ab_file)
    missing = [s for s in SAMPLES if s not in all_samples]
    if missing:
        sys.exit(f"samples not found in source table: {missing}")

    # Relative abundance: per-sample percentage of the KO counts shown here.
    shown = set(mat) & {k for k, v in ko2cl.items() if set(v) & KEEP}
    totals = {s: sum(mat[k][s] for k in shown) for s in SAMPLES}
    rel = {k: {s: mat[k][s] / totals[s] * 100 for s in SAMPLES} for k in shown}

    # Selection rule: up to MAX_PER_CLASS most abundant KOs per class, then cap.
    by_class = defaultdict(list)
    for ko in shown:
        for cl in ko2cl[ko]:
            if cl in KEEP:
                by_class[cl].append(ko)
    picked = set()
    for cl in KEEP:
        cand = sorted(set(by_class[cl]), key=lambda k: -sum(rel[k].values()))
        picked.update(cand[:MAX_PER_CLASS])
    mean_ab = {k: sum(rel[k].values()) / len(SAMPLES) for k in picked}
    if len(picked) > MAX_KO:
        picked = set(sorted(picked, key=lambda k: -mean_ab[k])[:MAX_KO])
    kos = sorted(picked)

    # ---- tables -------------------------------------------------------------
    with open(PROC / "ko_selected.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "KO",
                *[f"rel_abund_{s}" for s in SAMPLES],
                "mean_rel_abund",
                "classes",
                "n_pathways",
            ]
        )
        for ko in kos:
            cls = [c for c in ko2cl[ko] if c in KEEP]
            w.writerow(
                [
                    ko,
                    *[f"{rel[ko][s]:.6f}" for s in SAMPLES],
                    f"{mean_ab[ko]:.6f}",
                    ";".join(cls),
                    len(ko2path.get(ko, ())),
                ]
            )

    with open(PROC / "ko_links.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["KO", "class_code", "class_name", "connection_degree", "group"])
        for ko in kos:
            deg = len(ko2path.get(ko, ()))
            for cl in ko2cl[ko]:
                if cl in KEEP:
                    name = dict(CLASSES)[cl]
                    grp = "Metabolism" if cl in METABOLISM else "Other"
                    w.writerow([ko, cl, name, deg, grp])

    with open(PROC / "class_summary.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "class_code",
                "class_name",
                "n_links",
                "n_ko_total_in_kegg",
                *[f"sum_rel_abund_{s}" for s in SAMPLES],
                "group",
            ]
        )
        for code, name in CLASSES:
            members = [k for k in kos if code in ko2cl[k]]
            w.writerow(
                [
                    code,
                    name,
                    len(members),
                    len(by_class[code]),
                    *[f"{sum(rel[k][s] for k in members):.6f}" for s in SAMPLES],
                    "Metabolism" if code in METABOLISM else "Other",
                ]
            )

    prov = {
        "figure": "KEGG KO circos demo",
        "sources": [
            {
                "role": "KO -> pathway class hierarchy (level-1/2/3 + KO of ko00001)",
                "url": KEGG_URL,
                "file": str(kegg_file.relative_to(ROOT)),
                "sha256": sha256(kegg_file),
                "bytes": kegg_file.stat().st_size,
                "note": "KEGG REST snapshot; KEGG data are (c) Kanehisa Laboratories, academic use.",
            },
            {
                "role": "KEGG ortholog abundance profile (94 columns x 9437 KO rows)",
                "url": "https://figshare.com/articles/dataset/KEGG_KO_abundance_Mollisols_MR_SR_CT/29645648",
                "download": FIGSHARE_FILE,
                "doi": FIGSHARE_DOI,
                "license": "CC BY 4.0",
                "file": str(ab_file.relative_to(ROOT)),
                "sha256": sha256(ab_file),
                "bytes": ab_file.stat().st_size,
                "note": "Mollisols metagenomes; sample id = site letter + treatment letter + replicate.",
            },
        ],
        "derivation": {
            "samples_used": SAMPLES,
            "sample_meaning": "site L; A = chemical fertilizer only, B = manure, C = straw amendment",
            "n_samples_in_source": len(all_samples),
            "n_ko_in_source": len(mat),
            "n_ko_with_class_in_scope": len(shown),
            "n_ko_selected": len(kos),
            "selection_rule": f"top {MAX_PER_CLASS} KOs per class by mean relative abundance, cap {MAX_KO}, sorted by KO id",
            "relative_abundance": "per-sample percentage over the KO rows that belong to the 17 plotted classes",
            "simulated_values": False,
        },
    }
    (ROOT / "data" / "provenance.json").write_text(
        json.dumps(prov, indent=2), encoding="utf-8"
    )

    print(f"KOs in source          : {len(mat)}")
    print(f"KOs in plotted classes : {len(shown)}")
    print(f"KOs selected           : {len(kos)}")
    print(
        f"links (KO-class pairs) : {sum(1 for k in kos for c in ko2cl[k] if c in KEEP)}"
    )
    print(f"sample totals (counts) : {totals}")


if __name__ == "__main__":
    main()
