# Innovation Fair BD 2026 — 3 × 4 ft poster

A clinician-facing presentation of the decision support system: what a doctor sees and how
it helps them choose a dose. There is deliberately no mention of the inference engine, the
rule bases or the feedback loop.

| File | What it is |
|---|---|
| `FES-poster-3x4ft.pdf` | **Print this.** Vector, one page, exactly 914.4 × 1219.2 mm (3 ft wide × 4 ft tall). Sharp at any size. |
| `FES-poster-3x4ft-150dpi.png` | 5400 × 7200 px at 150 dpi, the usual resolution for large-format posters, for printers that require a raster. |
| `template.html` | Layout and copy. Edit this. |
| `build_poster.py` | Fills the template from a real run in the database. |
| `qr.swift` | Makes and reads QR codes with macOS's built-in Core Image. |
| `poster.html` | Generated — do not edit by hand, it is overwritten. |
| `proof.png` | Screen-size proof for review. |

## Before you print

- **Logo.** `poster/logo.png` (the University of Barishal logo) is embedded automatically;
  replace the file and rebuild to change it. Note its motto text is white, so it does not
  show against the white footer — only the emblem and the university name do.
- **GitHub Pages must be on**, or the QR code leads to a 404:
  repository **Settings → Pages → Branch `master`, folder `/docs` → Save**.
  Then open `https://moon1570.github.io/FES-tool/` once to confirm.

## Rebuilding

```bash
cd FESapp && ../.venv/bin/python ../poster/build_poster.py && cd ../poster
CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

# vector PDF at exactly 3 x 4 ft
"$CHROME" --headless=new --disable-gpu --no-pdf-header-footer \
  --virtual-time-budget=6000 --print-to-pdf=FES-poster-3x4ft.pdf "file://$PWD/poster.html"

# 150 dpi PNG (3456 CSS px x 1.5625 = 5400 px = 36 in)
"$CHROME" --headless=new --disable-gpu --hide-scrollbars --force-device-scale-factor=1.5625 \
  --virtual-time-budget=6000 --window-size=3456,4608 \
  --screenshot=FES-poster-3x4ft-150dpi.png "file://$PWD/poster.html"

# confirm the printed QR code still decodes to the landing page
swift qr.swift read FES-poster-3x4ft-150dpi.png
```

**If the PDF comes out as two pages**, the content has grown taller than 1219 mm. Reduce
the chart heights in `build_poster.py` or the section padding in `template.html`.

## Where the content comes from

Everything in the case study comes from the seeded run *Case 01 — standard*, read through the
same function that builds the website's home page (`app.views._case_study`), so the poster
and the site always show the same numbers:

- **10 billion → fewer than 1** tumour cells; without treatment the same tumour grows to
  about 18 billion (the model with every dose set to zero).
- **83 instead of 97** peak side-effect load, compared with the system's first suggestion.
- **4 of 9 doses lowered** before being given; 7% less drug in total.
- Lung, heart and kidney drug levels after each dose, with each organ's safe limit.

The QR code encodes `https://moon1570.github.io/FES-tool/` with error correction level Q
(still scans with roughly a quarter of it damaged), and the build is checked by reading the
code back from the rendered poster.

Typography for viewing from 2–3 m: title 44 mm, case-study figures 16.5 mm, body 9–10 mm.
Palette: blue `#1B4EE0` + coral `#E8543F` on white, `#16202B` text.
