"""Step 2: draw the three-layer circular figure with matplotlib only.

Layers, inside -> outside
  1. chord links      KO sector -> pathway class sector, line width = KEGG pathway count
  2. KO ring          three concentric heatmap rows, one per sample
  3. class arc        track, three grouped bars (summed relative abundance), label

The canvas is fitted to the drawn content, so the aspect ratio and the type sizes can be
changed without clipping labels or baking dead space into the figure.
"""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Circle, PathPatch, Rectangle, Wedge
from matplotlib.path import Path as MPath

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
OUT = ROOT / "output"

FIG_W, FIG_H = 11.2, 8.8  # ~1.27 : 1, chosen for the radial fan and for article display
TITLE_A = "KEGG ortholog abundance across"
TITLE_B = "three soil amendment treatments"
SUBTITLE = "Mollisols metagenomes   ·   48 KOs   ·   92 links   ·   17 KEGG classes"
NOTE_1 = ("bar height = summed relative abundance of the 48 KOs      "
          "link width = KEGG pathways per KO")
NOTE_2 = ("class arc = number of KOs      heatmap colour clipped at the "
          "5th / 95th percentile")
CREDIT = ("data: figshare 10.6084/m9.figshare.29645648 (CC BY 4.0)"
          "   ·   KEGG BRITE ko00001")

SAMPLE_NAMES = ["LA1", "LB1", "LC1"]
SAMPLE_TAG = ["control", "manure", "straw"]
SAMPLE_COLORS = ["#2F6E9E", "#D9A21B", "#8A63C4"]

GROUP_FILL = {"Metabolism": "#D6E4F1", "Other": "#F2CFDC"}
GROUP_EDGE = {"Metabolism": "#6F94B8", "Other": "#B96A82"}
LINK_COLOR = {"Metabolism": "#E8A227", "Other": "#5B90C6"}
INK = "#141414"
MUTED = "#6B6B6B"

CMAP = LinearSegmentedColormap.from_list(
    "abund", ["#2E9DA8", "#7FC6CC", "#F1EBED", "#F08CAF", "#D8326C"]
)

R_KO_IN = 0.560
R_LABEL = 0.594
R_HEAT = 0.742
ROW_H = 0.036
R_KO_OUT = R_HEAT + 3 * ROW_H
R_BAR_BASE = 0.660
BAR_MAX = 0.250
R_ARC = R_BAR_BASE + BAR_MAX + 0.012
R_PALE = R_BAR_BASE - 0.020
R_CLASS_LABEL = R_ARC + 0.045
TIER = 0.0             # single tier: thin type no longer needs the stagger

KO_SPAN, KO_START = 236.0, 298.0
CL_SPAN, CL_START = 112.0, -56.0
KO_GAP = 0.60
CLASS_GAP = 2.2

LEG_SCALE = 0.78          # bottom legend block: 1.0 = original size, 2/3 = compact
KO_FONT = 4.8
LX = -1.34             # left margin line
CLASS_FONT = 5.8          # one step above KO_FONT, same thin weight


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


def chord(ax, t0, r0, t1, r1, color, lw, alpha, gid=None):
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
            gid=gid,
        )
    )


def render(ax, kos, links, classes, ko_ang, cl_ang, norm, max_deg):
    ax.clear()
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.40, 1.95)
    ax.set_ylim(-1.35, 1.45)

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
                lw=0.9,
                zorder=0,
                gid=f'track-{c["class_code"]}',
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
                    lw=0.45,
                    zorder=3,
                    gid=f'bar-{c["class_code"]}-{i}',
                )
            )

    for ki, k in enumerate(kos):
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
                    lw=0.4,
                    zorder=3,
                    gid=f"cell-{ki}-{i}",
                )
            )

    ko_labels = [
        radial_text(
            ax,
            sum(ko_ang[k["KO"]]) / 2,
            R_LABEL,
            k["KO"],
            fontsize=KO_FONT,
            color="#2A2A2A",
            zorder=4,
            gid=f"kolab-{i}",
        )
        for i, k in enumerate(kos)
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
        lw = 0.4 + 1.9 * int(l["connection_degree"]) / max_deg
        chord(
            ax,
            end_frac[(i, False)],
            R_KO_IN,
            end_frac[(i, True)],
            R_BAR_BASE - 0.014,
            LINK_COLOR[l["group"]],
            lw,
            0.48,
            gid=f"chord-{i}",
        )

    for i, c in enumerate(classes):
        t0, t1 = cl_ang[c["class_code"]]
        tm = (t0 + t1) / 2
        r_lab = R_CLASS_LABEL + (TIER if i % 2 else 0.0)
        p0, p1 = xy(tm, R_ARC + 0.006), xy(tm, r_lab - 0.012)
        ax.plot(
            [p0[0], p1[0]],
            [p0[1], p1[1]],
            color=GROUP_EDGE[c["group"]],
            lw=0.9,
            zorder=2,
            solid_capstyle="butt",
        )
        radial_text(
            ax,
            tm,
            r_lab,
            c["class_name"],
            fontsize=CLASS_FONT,
            color="#2A2A2A",
            zorder=5,
            gid=f'clab-{c["class_code"]}',
        )

    # editorial header, top left
    ax.text(LX, 1.412, TITLE_A, fontsize=15.5, fontweight="bold", color=INK,
            ha="left", va="top", gid="head-title")
    ax.text(LX, 1.328, TITLE_B, fontsize=15.5, fontweight="bold", color=INK,
            ha="left", va="top")
    ax.text(LX, 1.178, SUBTITLE, fontsize=9.4, color="#3F3F3F", ha="left", va="top", gid="head-sub")
    ax.text(LX, 1.088, NOTE_1, fontsize=7.0, color=MUTED, ha="left", va="top", gid="head-note")
    ax.text(LX, 1.048, NOTE_2, fontsize=7.0, color=MUTED, ha="left", va="top")

    # legend strip, bottom — every dimension is multiplied by LEG_SCALE
    K = LEG_SCALE
    for i, (nm, tag) in enumerate(zip(SAMPLE_NAMES, SAMPLE_TAG)):
        y = -0.760 - i * 0.062
        ax.add_patch(Rectangle((LX, y), 0.040 * K, 0.040 * K,
                               facecolor=SAMPLE_COLORS[i], edgecolor="none", zorder=6, gid="leg"))
        ax.text(LX + 0.062, y + 0.020 * K, nm, fontsize=8.6 * K,
                fontweight="bold", va="center", color=INK)
        ax.text(LX + 0.230, y + 0.020 * K, tag, fontsize=8.6 * K,
                va="center", color="#3A3A3A", gid="leg")

    x0, y0, w, h = LX, -1.005, 0.52 * K, 0.026 * K
    for i in range(160):
        v = norm.inverse(0.02 + 0.96 * i / 159)
        ax.add_patch(Rectangle((x0 + w * i / 160, y0), w / 160 * 1.07 + 0.0006, h,
                               facecolor=CMAP(norm(v)), edgecolor="none", zorder=6, gid="cb"))
    ax.add_patch(Rectangle((x0, y0), w, h, facecolor="none",
                           edgecolor="#C9C9C9", lw=0.5, zorder=7, gid="cb"))
    for label, frac in ((f"{norm.vmin:.2f}%", 0.0), (f"{norm.vcenter:.2f}%", 0.5),
                        (f"{norm.vmax:.2f}%", 1.0)):
        px = x0 + w * frac
        ax.plot([px, px], [y0 - 0.014 * K, y0], color=INK, lw=0.6, zorder=7, gid="cbt")
        ax.text(px, y0 - 0.024 * K, label, fontsize=7.0 * K, ha="center",
                va="top", color="#2A2A2A")
    ax.text(x0, y0 + h + 0.014 * K, "Relative abundance per sample", fontsize=7.8 * K,
            fontweight="bold", va="bottom", color=INK, gid="cbt")

    sx, sy, sw = -0.70, y0, 0.44 * K
    ax.plot([sx, sx + sw], [sy, sy], color=INK, lw=0.8, zorder=6, gid="deg")
    for v in range(0, max_deg + 1, 5):
        px = sx + sw * v / max_deg
        ax.plot([px, px], [sy, sy + 0.020 * K], color=INK, lw=0.7, zorder=6, gid="deg")
        ax.text(px, sy - 0.022 * K, str(v), fontsize=7.0 * K, ha="center",
                va="top", color="#2A2A2A")
    ax.text(sx, sy + 0.034 * K, "Connection degree", fontsize=7.8 * K,
            fontweight="bold", va="bottom", color=INK, gid="deg")

    ax.text(LX, -1.085, CREDIT, fontsize=6.8 * K, color=MUTED, ha="left", va="top", gid="cred")

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

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    lims = None
    for _ in range(3):
        labels = render(ax, kos, links, classes, ko_ang, cl_ang, norm, max_deg)
        fig.canvas.draw()
        lims = content_box(ax, fig.canvas.get_renderer())
        pad = 0.012 * max(lims[1] - lims[0], lims[3] - lims[2])
        ax.set_xlim(lims[0] - pad, lims[1] + pad)
        ax.set_ylim(lims[2] - pad, lims[3] + pad)

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()

    worst = 0.0
    for t in labels:
        bb = t.get_window_extent(renderer=r).transformed(inv)
        worst = max(
            worst,
            max(math.hypot(cx, cy) for cx in (bb.x0, bb.x1) for cy in (bb.y0, bb.y1)),
        )
    print(
        f"KO label band: outermost corner r={worst:.3f}, heatmap inner r={R_HEAT:.3f}"
    )
    if worst > R_HEAT - 0.004:
        print("WARNING: KO labels touch the heatmap ring")

    for ext in ("png", "svg", "pdf"):
        fig.savefig(
            OUT / f"figure.{ext}",
            dpi=300,
            facecolor="white",
            bbox_inches="tight",
            pad_inches=0.10,
        )
    print(
        f"canvas x={tuple(round(v, 2) for v in ax.get_xlim())} "
        f"y={tuple(round(v, 2) for v in ax.get_ylim())}"
    )
    print(
        f"saved {OUT / 'figure.png'}  KOs={len(kos)} links={len(links)} "
        f"classes={len(classes)} max_degree={max_deg} range={vmin:.2f}-{vmax:.2f}%"
    )


if __name__ == "__main__":
    draw()
