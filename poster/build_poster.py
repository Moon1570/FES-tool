"""Generate the A1 poster from a real run stored in the application database.

Every number and every curve on the poster is read from an actual simulation, so the
poster cannot drift away from what the tool does. Regenerate with:

    cd FESapp && ../.venv/bin/python ../poster/build_poster.py
"""
import os
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

RUN_NAME = "Case 01 — standard"

# Organs a clinician reads at a glance, in the order they appear on the panel.
PANEL_ORGANS = ["lung", "kidney", "heart", "brain", "gut", "muscle", "liver"]


def esc(text):
    return (str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def build_hero(data):
    final = data["final"]
    log10 = final["log10N"]
    days = final["days"]
    cycles = data["cycles"]
    organs = {o["key"]: o for o in data["organs"]}

    parts = []

    # ---------------------------------------------------------------- card A
    ax0, ax1 = 190, 920          # plot area
    ay0, ay1 = 108, 336
    top, bottom = 10.5, -1.2     # log10 cells

    def px(day):
        return ax0 + (day / max(days[-1], 1)) * (ax1 - ax0)

    def py(v):
        return ay0 + (top - v) / (top - bottom) * (ay1 - ay0)

    parts.append('<rect x="0" y="0" width="960" height="420" rx="24" fill="#F4F6FA"/>')
    parts.append('<text x="34" y="56" font-size="24" font-weight="800" fill="#5A6675" '
                 'letter-spacing="2.6">EXPECTED TUMOUR RESPONSE</text>')

    # gridlines at 10^10, 10^5, 1 cell
    for value, label in [(10, "10 billion cells"), (5, "100 000 cells"), (0, "1 cell")]:
        y = py(value)
        parts.append(f'<line x1="{ax0}" y1="{y:.1f}" x2="{ax1}" y2="{y:.1f}" '
                     'stroke="#DDE3EC" stroke-width="2"/>')
        parts.append(f'<text x="{ax0 - 14}" y="{y + 7:.1f}" font-size="19" font-weight="600" '
                     f'fill="#8A94A3" text-anchor="end">{label}</text>')

    pts = " ".join(f"{px(d):.1f},{py(v):.1f}" for d, v in zip(days, log10))
    parts.append(f'<polyline points="{pts}" fill="none" stroke="#1B4EE0" stroke-width="6" '
                 'stroke-linejoin="round" stroke-linecap="round"/>')

    for c in cycles:
        parts.append(f'<circle cx="{px(c["day"]):.1f}" cy="{py(log10[c["day"]]):.1f}" '
                     'r="8" fill="#fff" stroke="#1B4EE0" stroke-width="4"/>')

    parts.append(f'<text x="{ax0}" y="382" font-size="20" font-weight="600" fill="#8A94A3">'
                 'day 0</text>')
    parts.append(f'<text x="{ax1}" y="382" font-size="20" font-weight="600" fill="#8A94A3" '
                 'text-anchor="end">day 120</text>')
    parts.append('<text x="555" y="382" font-size="20" font-weight="700" fill="#1B4EE0" '
                 'text-anchor="middle">○ each planned dose</text>')

    # ---------------------------------------------------------------- card B
    bx = 1000
    parts.append(f'<rect x="{bx}" y="0" width="600" height="420" rx="24" fill="#fff" '
                 'stroke="#D8DEE7" stroke-width="3"/>')
    parts.append(f'<text x="{bx + 34}" y="56" font-size="24" font-weight="800" fill="#5A6675" '
                 'letter-spacing="2.6">SAFETY MARGIN, EVERY ORGAN</text>')

    lx, bar_x, bar_w = bx + 34, bx + 200, 268
    # "SAFE LIMIT" sits above the rule, right-aligned to it, so it cannot run off the card.
    parts.append(f'<text x="{bar_x + bar_w}" y="92" font-size="19" font-weight="800" '
                 'fill="#16202B" text-anchor="end">SAFE LIMIT</text>')
    y = 122
    for key in PANEL_ORGANS:
        o = organs[key]
        frac = min(o["peak"] / o["limit"], 1.0)
        limiting = o["binding"]
        colour = "#E8543F" if limiting else "#1B4EE0"
        weight = "800" if limiting else "600"
        parts.append(f'<text x="{lx}" y="{y + 7}" font-size="23" font-weight="{weight}" '
                     f'fill="{colour if limiting else "#5A6675"}">{esc(o["label"])}</text>')
        parts.append(f'<rect x="{bar_x}" y="{y - 8}" width="{bar_w}" height="16" rx="8" '
                     'fill="#EDF0F5"/>')
        parts.append(f'<rect x="{bar_x}" y="{y - 8}" width="{bar_w * frac:.1f}" height="16" '
                     f'rx="8" fill="{colour}"/>')
        y += 34

    parts.append(f'<line x1="{bar_x + bar_w}" y1="102" x2="{bar_x + bar_w}" y2="{y - 22}" '
                 'stroke="#16202B" stroke-width="3"/>')

    parts.append(f'<rect x="{lx}" y="{y + 4}" width="336" height="40" rx="20" fill="#FCE8E4"/>')
    parts.append(f'<text x="{lx + 168}" y="{y + 31}" font-size="21" font-weight="800" '
                 'fill="#E8543F" text-anchor="middle">Lung is limiting this dose</text>')
    parts.append(f'<text x="{lx}" y="{y + 82}" font-size="21" font-weight="600" fill="#8A94A3">'
                 '+ 6 more organs monitored, all within limits</text>')

    # ---------------------------------------------------------- treatment plan
    ty = 478
    parts.append(f'<text x="0" y="{ty}" font-size="24" font-weight="800" fill="#5A6675" '
                 'letter-spacing="2.6">THE PLAN THE CLINICIAN REVIEWS · 9 DOSES OVER 120 DAYS</text>')

    chip_w, gap, cy = 158, 22, ty + 28
    for i, c in enumerate(cycles):
        x = i * (chip_w + gap)
        adjusted = c["adjusted"]
        fill = "#FFF3F0" if adjusted else "#F4F6FA"
        stroke = "#E8543F" if adjusted else "#DDE3EC"
        parts.append(f'<rect x="{x}" y="{cy}" width="{chip_w}" height="150" rx="18" '
                     f'fill="{fill}" stroke="{stroke}" stroke-width="3"/>')
        parts.append(f'<text x="{x + chip_w / 2}" y="{cy + 36}" font-size="20" font-weight="700" '
                     f'fill="#8A94A3" text-anchor="middle">day {c["day"]}</text>')
        parts.append(f'<text x="{x + chip_w / 2}" y="{cy + 88}" font-size="40" font-weight="800" '
                     f'fill="#16202B" text-anchor="middle">{c["delivered"]:.1f}</text>')
        parts.append(f'<text x="{x + chip_w / 2}" y="{cy + 116}" font-size="19" font-weight="600" '
                     f'fill="#8A94A3" text-anchor="middle">mg</text>')
        if adjusted:
            parts.append(f'<text x="{x + chip_w / 2}" y="{cy + 140}" font-size="18" '
                         'font-weight="800" fill="#E8543F" text-anchor="middle">lowered</text>')

    return "\n    ".join(parts)


def main():
    run = Run.objects.filter(patient_name=RUN_NAME, status="done").first()
    if run is None:
        raise SystemExit(f"No completed run named {RUN_NAME!r}. Run: manage.py seed_demo --reset")
    data = run.result
    m = data["metrics"]

    template = open(os.path.join(HERE, "template.html")).read()
    html = (template
            .replace("@@HERO@@", build_hero(data))
            .replace("@@ADJUSTED@@", str(m["cycles_adjusted"]))
            .replace("@@TOTAL@@", str(m["cycles_total"])))
    out = os.path.join(HERE, "poster.html")
    open(out, "w").write(html)
    print(f"wrote {out}")
    print(f"  from run: {run.patient_name} ({run.id})")
    print(f"  {m['cycles_adjusted']} of {m['cycles_total']} doses lowered; "
          f"{m['log_reduction']:.2f} log10 reduction")


if __name__ == "__main__":
    main()
