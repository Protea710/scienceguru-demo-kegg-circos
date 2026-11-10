"""Step 2: draw the three-layer circular figure with matplotlib only.

Layers, inside -> outside
  1. chord links      KO sector -> pathway class sector, line width = KEGG pathway count
  2. KO ring          three concentric heatmap rows, one per sample
  3. class arc        pale track, three grouped bars (summed relative abundance), label

The canvas is fitted to the drawn content: artists are laid out once, their bounding box
is measured, and the axis limits are set from it, so nothing is clipped and no dead
space is baked into the figure.
"""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import PathPatch, Rectangle, Wedge
from matplotlib.path import Path as MPath

ROOT = Path(__file__).resolve().parents[2]   # repo root: this file lives in archive/v1/
PROC = ROOT / "data" / "processed"
OUT = Path(__file__).resolve().parent    # keep v1 output next to this script

SAMPLE_NAMES = ["LA1", "LB1", "LC1"]
SAMPLE_TAG = ["control", "manure", "straw"]
# outside the heatmap's teal-pink range, so colour never carries two meanings
SAMPLE_COLORS = ["#2F6E9E", "#D9A21B", "#8A63C4"]
GROUP_FILL = {"Metabolism": "#E6EFF7", "Other": "#FAE9EE"}
GROUP_EDGE = {"Metabolism": "#9DB9D0", "Other": "#E7A9BA"}
LINK_COLOR = {"Metabolism": "#E8A227", "Other": "#5B90C6"}
INK = "#1b1b1b"
MUTED = "#6d6d6d"

CMAP = LinearSegmentedColormap.from_list(
    "abund", ["#2E9DA8", "#93CDD3", "#FAF4F6", "#F49CBB", "#DC3E76"]
)

R_KO_IN = 0.560
R_LABEL = 0.606
R_HEAT = 0.724
ROW_H = 0.036
R_KO_OUT = R_HEAT + 3 * ROW_H
R_BAR_BASE = 0.660
BAR_MAX = 0.250
R_ARC = R_BAR_BASE + BAR_MAX + 0.012
R_PALE = R_BAR_BASE - 0.020
R_CLASS_LABEL = R_ARC + 0.045
TIER = 0.105

KO_SPAN, KO_START = 236.0, 298.0
CL_SPAN, CL_START = 112.0, -56.0
KO_GAP = 0.55
CLASS_GAP = 1.8

KO_FONT = 4.6
CLASS_FONT = 7.8
LX = -1.15  # legend column


def xy(t, r):
    a = math.radians(t)
    return r * math.cos(a), r * math.sin(a)


def load():
    kos = list(
        csv.DictReader(open(PROC / "ko_selected.tsv", encoding="utf-8"), delimiter="\t")
    )
    links = list(
        csv.DictReader(open(PROC / "ko_links.tsv", encoding="utf-8"), delimiter="\t")
    )
    classes = list(
        csv.DictReader(
            open(PROC / "class_summary.tsv", encoding="utf-8"), delimiter="\t"
        )
    )
    return kos, links, classes


def build_angles(kos, classes):
    n_ko = len(kos)
    w_ko = (KO_SPAN - KO_GAP * n_ko) / n_ko
    ko_ang, t = {}, KO_START
    for k in kos:
        ko_ang[k["KO"]] = (t - w_ko, t)
        t -= w_ko + KO_GAP

    weights = [max(int(c["n_links"]), 1) for c in classes]
    usable = CL_SPAN - CLASS_GAP * len(classes)
    cl_ang, t = {}, CL_START
    for c, w in zip(classes, weights):
        span = usable * w / sum(weights)
        cl_ang[c["class_code"]] = (t, t + span)
        t += span + CLASS_GAP
    return ko_ang, cl_ang


def radial_text(ax, t, r, s, **kw):
    x, y = xy(t, r)
    rot, ha = t, "left"
    if 90 < (t % 360) < 270:
        rot, ha = t + 180, "right"
    kw.setdefault("va", "center")
    return ax.text(x, y, s, rotation=rot, rotation_mode="anchor", ha=ha, **kw)


def chord(ax, t0, r0, t1, r1, color, lw, alpha):
    p0, p1 = xy(t0, r0), xy(t1, r1)
    c0, c1 = xy(t0, r0 * 0.34), xy(t1, r1 * 0.34)
    path = MPath(
        [p0, c0, c1, p1], [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
    )
    ax.add_patch(
        PathPatch(
            path,
            facecolor="none",
            edgecolor=color,
            lw=lw,
            alpha=alpha,
            capstyle="round",
            zorder=1,
        )
    )


def render(ax, kos, links, classes, ko_ang, cl_ang, norm, max_deg):
    ax.clear()
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.30, 1.90)
    ax.set_ylim(-1.30, 1.15)

    bar_max = max(float(c[f"sum_rel_abund_{s}"]) for c in classes for s in SAMPLE_NAMES)
    for c in classes:
        t0, t1 = cl_ang[c["class_code"]]
        ax.add_patch(
            Wedge(
                (0, 0),
                R_ARC,
                t0,
                t1,
                width=R_ARC - R_PALE,
                facecolor=GROUP_FILL[c["group"]],
                edgecolor=GROUP_EDGE[c["group"]],
                lw=0.5,
                zorder=0,
            )
        )

    for c in classes:
        t0, t1 = cl_ang[c["class_code"]]
        span = t1 - t0
        bw = span / 3.0 * 0.90
        for i, s in enumerate(SAMPLE_NAMES):
            v = float(c[f"sum_rel_abund_{s}"])
            a0 = t0 + span * (0.05 + i / 3.0)
            h = BAR_MAX * v / bar_max
            ax.add_patch(
                Wedge(
                    (0, 0),
                    R_BAR_BASE + h,
                    a0,
                    a0 + bw,
                    width=h,
                    facecolor=SAMPLE_COLORS[i],
                    edgecolor="white",
                    lw=0.3,
                    zorder=3,
                )
            )

    for k in kos:
        t0, t1 = ko_ang[k["KO"]]
        for i, s in enumerate(SAMPLE_NAMES):
            r0 = R_HEAT + (2 - i) * ROW_H
            ax.add_patch(
                Wedge(
                    (0, 0),
                    r0 + ROW_H * 0.97,
                    t0,
                    t1,
                    width=ROW_H * 0.94,
                    facecolor=CMAP(norm(float(k[f"rel_abund_{s}"]))),
                    edgecolor="white",
                    lw=0.35,
                    zorder=3,
                )
            )

    ko_labels = [
        radial_text(
            ax,
            sum(ko_ang[k["KO"]]) / 2,
            R_LABEL,
            k["KO"],
            fontsize=KO_FONT,
            color=INK,
            zorder=4,
        )
        for k in kos
    ]

    by_class, by_ko = {}, {}
    for i, l in enumerate(links):
        by_class.setdefault(l["class_code"], []).append(i)
        by_ko.setdefault(l["KO"], []).append(i)
    end_frac = {}
    for store, ang, side in ((by_class, cl_ang, True), (by_ko, ko_ang, False)):
        for key, idxs in store.items():
            t0, t1 = ang[key]
            for j, idx in enumerate(idxs):
                end_frac[(idx, side)] = t0 + (t1 - t0) * (
                    0.12 + 0.76 * (j + 0.5) / len(idxs)
                )

    for i, l in enumerate(links):
        lw = 0.35 + 1.9 * int(l["connection_degree"]) / max_deg
        chord(
            ax,
            end_frac[(i, False)],
            R_KO_IN,
            end_frac[(i, True)],
            R_BAR_BASE - 0.014,
            LINK_COLOR[l["group"]],
            lw,
            0.46,
        )

    for i, c in enumerate(classes):
        t0, t1 = cl_ang[c["class_code"]]
        tm = (t0 + t1) / 2
        r_lab = R_CLASS_LABEL + (TIER if i % 2 else 0.0)
        p0, p1 = xy(tm, R_ARC + 0.006), xy(tm, r_lab - 0.010)
        ax.plot(
            [p0[0], p1[0]],
            [p0[1], p1[1]],
            color=GROUP_EDGE[c["group"]],
            lw=0.7,
            zorder=2,
        )
        radial_text(
            ax,
            tm,
            r_lab,
            c["class_name"],
            fontsize=CLASS_FONT,
            color="#141414",
            zorder=5,
        )

    # legends: one header line on top, one strip at the bottom, both clear of the ring
    ax.text(-1.16, 0.995,
            f"{len(kos)} KEGG orthologs   ·   {len(links)} KO–class links   ·   "
            f"{len(classes)} KEGG classes",
            fontsize=8.6, ha="left", va="bottom", color=INK)
    ax.text(-1.16, 0.928, "bar height = summed relative abundance of the 48 KOs",
            fontsize=6.8, ha="left", va="bottom", color=MUTED)
    ax.text(-1.16, 0.888, "link width = KEGG pathways per KO  ·  class arc = number of KOs",
            fontsize=6.8, ha="left", va="bottom", color=MUTED)

    for i, (nm, tag) in enumerate(zip(SAMPLE_NAMES, SAMPLE_TAG)):
        x = -1.16 + i * 0.62
        ax.add_patch(Rectangle((x, -0.960), 0.034, 0.034,
                               facecolor=SAMPLE_COLORS[i], edgecolor="none", zorder=6))
        ax.text(x + 0.052, -0.943, f"{nm}   {tag}", fontsize=7.8, va="center", color=INK)

    x0, y0, w, h = -1.16, -1.145, 0.50, 0.026
    for i in range(160):
        v = norm.inverse(0.02 + 0.96 * i / 159)
        ax.add_patch(Rectangle((x0 + w * i / 160, y0), w / 160 * 1.07, h,
                               facecolor=CMAP(norm(v)), edgecolor="none", zorder=6))
    for label, frac in ((f"{norm.vmin:.2f}%", 0.0), (f"{norm.vcenter:.2f}%", 0.5),
                        (f"{norm.vmax:.2f}%", 1.0)):
        px = x0 + w * frac
        ax.plot([px, px], [y0 - 0.015, y0], color=INK, lw=0.7, zorder=7)
        ax.text(px, y0 - 0.026, label, fontsize=6.8, ha="center", va="top")
    ax.text(x0, y0 + h + 0.012, "Relative abundance per sample", fontsize=7.2, va="bottom")
    ax.text(x0, y0 - 0.072, "colour clipped at the 5th / 95th percentile",
            fontsize=6.2, va="top", color=MUTED)

    sx, sy, sw = 0.02, -1.145, 0.40
    ax.plot([sx, sx + sw], [sy, sy], color=INK, lw=0.9, zorder=6)
    for v in range(0, max_deg + 1, 5):
        px = sx + sw * v / max_deg
        ax.plot([px, px], [sy, sy + 0.020], color=INK, lw=0.8, zorder=6)
        ax.text(px, sy - 0.022, str(v), fontsize=6.6, ha="center", va="top")
    ax.text(sx, sy + 0.032, "Connection degree", fontsize=7.0, va="bottom")
    ax.text(sx, sy - 0.072, "KEGG pathways per ortholog", fontsize=6.2, va="top", color=MUTED)

    return ko_labels


def content_box(ax, renderer):
    inv = ax.transData.inverted()
    xs, ys = [-R_KO_OUT, R_KO_OUT], [-R_KO_OUT, R_KO_OUT]
    for t in ax.texts:
        bb = t.get_window_extent(renderer=renderer).transformed(inv)
        xs += [bb.x0, bb.x1]
        ys += [bb.y0, bb.y1]
    return min(xs), max(xs), min(ys), max(ys)


def draw():
    kos, links, classes = load()
    OUT.mkdir(exist_ok=True)
    ko_ang, cl_ang = build_angles(kos, classes)

    vals = sorted(float(k[f"rel_abund_{s}"]) for k in kos for s in SAMPLE_NAMES)
    n = len(vals)
    vmin, vmax, vmid = vals[int(n * 0.05)], vals[int(n * 0.95) - 1], vals[n // 2]
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vmid, vmax=vmax)
    max_deg = max(int(l["connection_degree"]) for l in links)

    fig, ax = plt.subplots(figsize=(10.6, 9.4))
    lims = None
    for _ in range(3):
        labels = render(ax, kos, links, classes, ko_ang, cl_ang, norm, max_deg)
        fig.canvas.draw()
        lims = content_box(ax, fig.canvas.get_renderer())
        pad = 0.015 * max(lims[1] - lims[0], lims[3] - lims[2])
        ax.set_xlim(lims[0] - pad, lims[1] + pad)
        ax.set_ylim(lims[2] - pad, lims[3] + pad)

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()

    worst = 0.0
    for t in labels:
        bb = t.get_window_extent(renderer=r).transformed(inv)
        worst = max(worst, max(math.hypot(cx, cy)
                               for cx in (bb.x0, bb.x1) for cy in (bb.y0, bb.y1)))
    print(f"KO label band: outermost corner r={worst:.3f}, heatmap inner r={R_HEAT:.3f}")
    if worst > R_HEAT - 0.004:
        print("WARNING: KO labels touch the heatmap ring")

    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"figure.{ext}", dpi=300, facecolor="white",
                    bbox_inches="tight", pad_inches=0.08)
    print(f"canvas limits x={tuple(round(v,2) for v in ax.get_xlim())} "
          f"y={tuple(round(v,2) for v in ax.get_ylim())}")
    print(f"saved {OUT/'figure.png'}  KOs={len(kos)} links={len(links)} "
          f"classes={len(classes)} max_degree={max_deg} range={vmin:.2f}-{vmax:.2f}%")



if __name__ == "__main__":
    draw()
