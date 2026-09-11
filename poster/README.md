# Innovation Fair BD 2026 — A1 poster

A clinician-facing presentation of the decision support system. It shows what a doctor
sees and how the tool helps them choose a dose — not how the software works. There is
deliberately no mention of the inference engine, the rule bases or the feedback loop.

| File | What it is |
|---|---|
| `FES-poster-A1.pdf` | **Print this.** Vector, one page, exactly 594 × 841 mm. Sharper than the raster at any size. |
| `FES-poster-A1-300dpi.png` | 7016 × 9934 px at 300 dpi, for printers that require a raster. |
| `template.html` | The layout and copy. Edit this. |
| `build_poster.py` | Fills the template from a real run in the database. |
| `poster.html` | Generated — do not edit by hand, it is overwritten. |
| `proof.png` | Screen-size proof for review. |

## Before you print — three things to fill in

All three are in the footer of `template.html`:

1. **Author names and institution** — replace the placeholder in
   `<div class="names">Author names — Institution / Department</div>`, then delete the
   `border` and `color` rules on `.team .names` so it stops looking like a form field.
2. **Institution logo** — drop your logo into the first `.slot` box.
3. **QR code** — generate one pointing at your demo or repository, into the second `.slot`.

## Regenerating

The hero is drawn from an actual stored run, so the poster cannot drift from what the
tool really produces. If the seeded data is missing, create it first with
`manage.py seed_demo --reset`.

```bash
cd FESapp && ../.venv/bin/python ../poster/build_poster.py && cd ../poster

CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

# vector PDF at exact A1
"$CHROME" --headless --disable-gpu --no-pdf-header-footer \
  --virtual-time-budget=7000 --print-to-pdf=FES-poster-A1.pdf "file://$PWD/poster.html"

# 300 dpi PNG (2245 CSS px x 3.125 = 7016 px = 594 mm at 300 dpi)
"$CHROME" --headless --disable-gpu --hide-scrollbars \
  --force-device-scale-factor=3.125 --virtual-time-budget=7000 \
  --window-size=2245,3179 --screenshot=FES-poster-A1-300dpi.png "file://$PWD/poster.html"
```

**If the PDF comes out as two pages**, the content has grown taller than 841 mm. Trim the
vertical padding on `.results`, `.process` or `.cols` until `document.documentElement.scrollHeight`
equals the body height (3179 px).

## Where the content comes from

Read from the seeded run *Case 01 — standard*:

- The response curve is that run's actual 121-day tumour trajectory, with a marker on
  each of the nine dose days.
- The organ panel shows seven organs' real peak concentrations against their limits, with
  the lung flagged because it is the organ that actually limited the dose in that run.
- The nine dose chips are the delivered doses; the four marked "lowered" are the four the
  safety check reduced.
- **13** organs checked, **4 of 9** doses lowered, **28** alternative schedules — all from
  the application, not estimates.

Typography is set for 2 m legibility: headline 28.5 mm, KPI numerals 60 mm, body 8 mm.
Palette is blue `#1B4EE0` + coral `#E8543F` on white with `#16202B` text.
