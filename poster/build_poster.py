"""Generate the 3 x 4 ft poster from a real run stored in the application database.

Every number and curve on the poster is read from an actual simulation, through the
same function that builds the home page's case study, so the poster cannot drift from
what the site shows. The QR code is made with macOS's built-in Core Image (qr.swift).

    cd FESapp && ../.venv/bin/python ../poster/build_poster.py
"""
import base64
import glob
import math
import os
import subprocess
import sys
import warnings

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.join(os.path.dirname(HERE), "FESapp")
sys.path.insert(0, APP)
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "FESapp.settings")

import django  # noqa: E402

django.setup()

from app.models import Run  # noqa: E402
from app.views import _case_study  # noqa: E402

RUN_NAME = "Case 01 — standard"
QR_URL = "https://moon1570.github.io/FES-tool/"
QR_TEXT = "moon1570.github.io/FES-tool"

BLUE, CORAL, AMBER, INK = "#1B4EE0", "#E8543F", "#B26A00", "#16202B"
GREY, MUTED, GRID = "#5A6675", "#8A94A3", "#DDE3EC"
DOSE_COLOURS = ["#1f6feb", "#2b8ae0", "#22a0c8", "#1fae9e", "#3fae63",
                "#8aa63a", "#c09428", "#cf6f24", "#c0392b", "#a3306e"]


def esc(text):
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def cells_in_words(log10):
    if log10 is None:
        return "—"
    if log10 < 0:
        return "fewer than 1"
    n = 10 ** log10
    for unit, word in ((1e9, "billion"), (1e6, "million"), (1e3, "thousand")):
        if n >= unit:
            x = n / unit
            return (f"{round(x)}" if x >= 10 else f"{x:.1f}".rstrip("0").rstrip(".")) + " " + word
    return str(round(n))


def nice_ticks(top, count=3):
    """Round tick values from 0 to at most `top`."""
    raw = top / count
    mag = 10 ** math.floor(math.log10(raw))
    step = next(m * mag for m in (1, 2, 5, 10) if m * mag >= raw)
    return [i * step for i in range(int(top // step) + 1)]


class Chart:
    """A minimal SVG line/bar chart in viewBox units, styled for print."""

    def __init__(self, w, h, box, xr, yr, fs=26):
        self.w, self.h, self.fs = w, h, fs
        self.l, self.r, self.t, self.b = box
        self.x0, self.x1 = xr
        self.y0, self.y1 = yr
        self.parts = []

    def px(self, x):
        return self.l + (x - self.x0) / (self.x1 - self.x0) * (self.w - self.l - self.r)

    def py(self, y):
        y = min(max(y, self.y0), self.y1)
        return self.h - self.b - (y - self.y0) / (self.y1 - self.y0) * (self.h - self.t - self.b)

    def axes(self, yticks, xticks):
        for v, label in yticks:
            y = self.py(v)
            self.parts.append(f'<line x1="{self.l}" y1="{y:.1f}" x2="{self.w - self.r}" y2="{y:.1f}" '
                              f'stroke="{GRID}" stroke-width="2"/>')
            self.parts.append(f'<text x="{self.l - 10}" y="{y + self.fs * 0.34:.1f}" font-size="{self.fs}" '
                              f'fill="{MUTED}" text-anchor="end">{esc(label)}</text>')
        base = self.h - self.b
        for v, label in xticks:
            x = self.px(v)
            anchor = "start" if v == self.x0 else ("end" if v == self.x1 else "middle")
            self.parts.append(f'<text x="{x:.1f}" y="{base + self.fs * 1.25:.1f}" font-size="{self.fs}" '
                              f'fill="{MUTED}" text-anchor="{anchor}">{esc(label)}</text>')

    def hline(self, y, colour, width=4, dash="14 10"):
        self.parts.append(f'<line x1="{self.l}" y1="{self.py(y):.1f}" x2="{self.w - self.r}" '
                          f'y2="{self.py(y):.1f}" stroke="{colour}" stroke-width="{width}" '
                          f'stroke-dasharray="{dash}"/>')

    def line(self, xs, ys, colour, width, dash=None, fill=None):
        pts = [(self.px(x), self.py(y)) for x, y in zip(xs, ys) if y is not None]
        if not pts:
            return
        d = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        if fill:
            base = self.py(self.y0)
            self.parts.append(f'<polygon points="{pts[0][0]:.1f},{base:.1f} {d} {pts[-1][0]:.1f},{base:.1f}" '
                              f'fill="{fill}"/>')
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        self.parts.append(f'<polyline points="{d}" fill="none" stroke="{colour}" stroke-width="{width}" '
                          f'stroke-linejoin="round" stroke-linecap="round"{dash_attr}/>')

    def bars(self, xs, ys, colour, bar_w):
        base = self.py(self.y0)
        half = (self.px(bar_w) - self.px(0)) / 2
        for x, y in zip(xs, ys):
            if y and y > 0:
                top = self.py(y)
                self.parts.append(f'<rect x="{self.px(x) - half:.1f}" y="{top:.1f}" width="{2 * half:.1f}" '
                                  f'height="{base - top:.1f}" fill="{colour}"/>')

    def svg(self, label):
        return (f'<svg viewBox="0 0 {self.w} {self.h}" role="img" aria-label="{esc(label)}" '
                f'font-family="Helvetica Neue, Helvetica, Arial, sans-serif">' + "".join(self.parts) + "</svg>")


def tumour_chart(c):
    days = c["days"]
    ch = Chart(800, 520, (170, 22, 18, 52), (days[0], days[-1]), (-1.5, 11.2))
    ch.axes([(10, "10 billion"), (6, "1 million"), (2, "100"), (0, "1")],
            [(0, "day 0"), (60, "day 60"), (days[-1], f"day {days[-1]}")])
    ch.line(days, c["untreated"], CORAL, 6, dash="16 11")
    ch.line(days, c["treated"], BLUE, 8)
    return ch.svg("Tumour cells over the course: no treatment against the recommended plan")


def toxicity_chart(c):
    days, limit = c["days"], c["metrics"]["toxicity_limit"]
    ch = Chart(800, 520, (82, 22, 18, 52), (days[0], days[-1]), (0, limit * 1.12))
    ch.axes([(0, "0"), (50, "50"), (100, "100")],
            [(0, "day 0"), (60, "day 60"), (days[-1], f"day {days[-1]}")])
    ch.hline(limit, CORAL, 4)
    ch.line(days, c["tox_planned"], "#7A8796", 4, dash="3 9")
    ch.line(days, c["tox_final"], AMBER, 6, fill="rgba(178,106,0,0.12)")
    return ch.svg("Side-effect load over the course: first suggestion against the recommended plan")


def dose_chart(c):
    days = c["days"]
    top = max(v for v in c["dose_planned"] if v is not None)
    ch = Chart(800, 520, (82, 22, 18, 52), (days[0] - 2, days[-1] + 2), (0, top * 1.12))
    ch.axes([(v, f"{v:g}") for v in nice_ticks(top * 1.05)],
            [(days[0] - 2, "day 0"), (60, "day 60"), (days[-1] + 2, f"day {days[-1]}")])
    ch.bars(days, c["dose_planned"], "rgba(122,135,150,0.40)", 3.0)
    ch.bars(days, c["dose_final"], BLUE, 3.0)
    return ch.svg("Dose given on each day: first suggestion against the recommended plan")


def organ_chart(o, hours):
    top = max(o["limit"], max(max(s) for s in o["series"])) * 1.12
    ch = Chart(800, 385, (64, 34, 14, 50), (0, hours[-1]), (0, top))
    ch.axes([(v, f"{v:g}") for v in nice_ticks(top * 0.98)],
            [(12, "12"), (24, "24"), (36, "36"), (hours[-1], "48 h")])
    ch.hline(o["limit"], CORAL, 4.5)
    for i, series in enumerate(o["series"]):
        ch.line(hours, series, DOSE_COLOURS[i % len(DOSE_COLOURS)], 4)
    return ch.svg(f"Drug level in the {o['label'].lower()} over 48 hours after each dose")


def qr_svg(url):
    out = subprocess.run(["swift", os.path.join(HERE, "qr.swift"), "gen", url],
                         capture_output=True, text=True, cwd=HERE)
    rows = [r for r in out.stdout.split() if r]
    if out.returncode != 0 or not rows:
        raise SystemExit(f"QR generation failed: {out.stderr.strip()}")
    n, q = len(rows), 4                       # 4-module quiet zone, as the QR spec requires
    size = n + 2 * q
    path = "".join(f"M{x + q},{y + q}h1v1h-1z"
                   for y, row in enumerate(rows) for x, ch in enumerate(row) if ch == "1")
    return (f'<svg viewBox="0 0 {size} {size}" shape-rendering="crispEdges" role="img" '
            f'aria-label="QR code for {esc(url)}"><rect width="{size}" height="{size}" fill="#fff"/>'
            f'<path d="{path}" fill="#000"/></svg>')


def logo_html():
    """Embed poster/logo.(svg|png|jpg) if present; otherwise leave a labelled space."""
    for path in sorted(glob.glob(os.path.join(HERE, "logo.*"))):
        ext = path.rsplit(".", 1)[-1].lower()
        mime = {"svg": "image/svg+xml", "png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg"}.get(ext)
        if mime:
            data = base64.b64encode(open(path, "rb").read()).decode()
            return f'<img src="data:{mime};base64,{data}" alt="University of Barishal logo">'
    return '<div class="slot">University<br>of Barishal<br>logo</div>'


def main():
    run = Run.objects.filter(patient_name=RUN_NAME, status="done").first()
    if run is None or not run.result:
        raise SystemExit(f"No completed run named {RUN_NAME!r}. Run: manage.py seed_demo --reset")
    c = _case_study(run)
    m = c["metrics"]
    burden = run.result["patient"].get("burden", {})
    less = 100 * (1 - m["total_delivered"] / m["total_planned"]) if m["total_planned"] else 0

    organs = []
    for o in c["organs"]:
        chip = ('<span class="chip cut">Limited the dose</span>' if o["binding"]
                else '<span class="chip ok">Within its limit</span>')
        organs.append(
            f'<div class="organ{" limiting" if o["binding"] else ""}">'
            f'<div class="organ-head"><strong>{esc(o["label"])}</strong>{chip}</div>'
            f'<div class="organ-sub">Highest level {o["peak"]:.1f} · safe limit {o["limit"]:.1f}</div>'
            f'{organ_chart(o, c["hours"])}</div>')
    dose_key = "".join(f'<span><i style="background:{DOSE_COLOURS[i % len(DOSE_COLOURS)]}"></i>Dose {n}</span>'
                       for i, n in enumerate(c["doses"]))

    fill = {
        "@@CASE_HEADING@@": esc(f"A {run.weight_kg:.0f} kg adult with {run.burden_label.lower()}, "
                                f"treated for {run.horizon_days} days"),
        "@@CASE_SUB@@": esc(f"About {cells_in_words(m['log10_n_start'])} tumour cells at the start "
                            f"({burden.get('mass', '')} of tumour), a dose every {run.interval_days} days, "
                            f"and {run.preset_label.lower()} organ limits. Every figure here is the "
                            f"system's own output for this case."),
        "@@P1_BIG@@": esc(f"{cells_in_words(m['log10_n_start'])} → "
                          + ("fewer than 1" if m["log10_n_end"] < 0 else cells_in_words(m["log10_n_end"]))),
        "@@P1_CAP@@": esc(f"tumour cells after {c['days'][-1]} days. Without treatment, the same tumour "
                          f"grows to about {cells_in_words(c['untreated'][-1])}."),
        "@@P1_SVG@@": tumour_chart(c),
        "@@P2_BIG@@": esc(f"{round(m['peak_toxicity'])} instead of {round(c['peak_tox_planned'])}"),
        "@@P2_CAP@@": esc("peak side-effect load with the recommended plan, against the system's "
                          "first suggestion."),
        "@@P2_SVG@@": toxicity_chart(c),
        "@@P3_BIG@@": esc(f"{m['cycles_adjusted']} of {m['cycles_total']} doses lowered"),
        "@@P3_CAP@@": esc("before being given, so every organ stays within its safe limit"
                          + (f" — {less:.0f}% less drug in total." if less else ".")),
        "@@P3_SVG@@": dose_chart(c),
        "@@DOSE_KEY@@": dose_key,
        "@@ORGANS@@": "".join(organs),
        "@@LOGO@@": logo_html(),
        "@@QR_SVG@@": qr_svg(QR_URL),
        "@@QR_TEXT@@": esc(QR_TEXT),
    }
    html = open(os.path.join(HERE, "template.html")).read()
    for key, value in fill.items():
        if key not in html:
            raise SystemExit(f"template is missing {key}")
        html = html.replace(key, value)
    out = os.path.join(HERE, "poster.html")
    open(out, "w").write(html)
    print(f"wrote {out}")
    print(f"  from run: {run.patient_name} ({run.id})")
    print(f"  {fill['@@P1_BIG@@']} | {fill['@@P2_BIG@@']} | {fill['@@P3_BIG@@']}")
    print(f"  organs: {', '.join(o['label'] for o in c['organs'])} | QR -> {QR_URL}")


if __name__ == "__main__":
    main()
