# FES-tool

**Decision Support System for Cancer Chemotherapy** — dose scheduling by fuzzy inference
under physiologically-based organ safety constraints.

A two-stage **Fuzzy Expert System** (FES-1 → FES-2) proposes a dose for each treatment
cycle from tumour burden, cumulative toxicity and body surface area. Every proposed dose
is then pushed through a **35-compartment PBPK model** and checked against per-organ peak
exposure ceilings; any dose that breaches a ceiling is reduced by 5% and re-checked until
it passes. The tumour model is then re-simulated with the corrected schedule.

> Simulation output, for research and evaluation. Not a medical device; not for clinical
> decision-making.

Python port of the MATLAB work behind
*Healthcare Analytics* (2023), [S2772442523000060](https://www.sciencedirect.com/science/article/pii/S2772442523000060).

## Running it

```bash
python3.11 -m venv .venv
.venv/bin/pip install -r FESapp/requirements.txt
cd FESapp
../.venv/bin/python manage.py migrate
../.venv/bin/python manage.py seed_demo --reset    # optional: 4 presentable runs
../.venv/bin/python manage.py runserver
```

Open <http://127.0.0.1:8000>. Everything runs locally — the Plotly bundle is vendored in
`static/js/`, so the app works with no network connection.

## Deploying

To put the app online for free, with no cold start, follow [DEPLOY.md](DEPLOY.md)
(PythonAnywhere for the app, GitHub Pages for the landing page the QR code points to).

## Layout

```
FESapp/
├── app/
│   ├── engine/            the science, importable and free of Django
│   │   ├── params.py      tumour + PBPK constants
│   │   ├── fuzzy.py       FES-1 / FES-2 (rule graphs built once at import)
│   │   ├── tumor.py       Martin's 3-state tumour/PK/toxicity model
│   │   ├── pbpk.py        35-state PBPK model; STATE_NAMES is the ordering authority
│   │   ├── safety.py      per-organ ceilings, limit presets, structured report
│   │   ├── burden.py      initial tumour burden presets and clinical descriptors
│   │   ├── sweep.py       regimen grid search and trade-off frontier
│   │   ├── _solver.py     serialised access to LSODA (it is not reentrant)
│   │   ├── runner.py      orchestration + the 5% dose-reduction repair loop
│   │   └── serialize.py   RunResult -> JSON payload for the UI
│   ├── executor.py        runs simulations on a worker thread, tracks progress
│   ├── models.py          Run: inputs, status and stored result
│   └── views.py           form handling, polling, dashboard, CSV export
├── templates/             base, run_setup, result, explore, sweep, history, compare
├── static/                app.css, charts.js, vendored plotly-basic
├── tools/regress.py       regression suite (see below)
└── tests/baseline_original.json   recorded output of the pre-rewrite code
```

## Verification

```bash
cd FESapp && ../.venv/bin/python tools/regress.py
```

The suite pins the rewrite against `tests/baseline_original.json`, a recording of the
original `views.calc`. In legacy mode (`log_tumor=False`, `legacy_dose_gate=True`,
50.0 ceilings) the new engine reproduces the original's dose schedule, BSA and N(84)
to within 1e-9 — while running about 350x faster.

## What changed, and why

### Defects fixed

| Defect | Effect |
|---|---|
| `dose` was a module global, never reset, and memoised per day | The second patient simulated in a process silently reused the first patient's schedule. Weight had no effect on the result at all. |
| The dose-vs-day chart plotted a hardcoded 121-element literal | Every user saw the same canned plot regardless of input. |
| `feedback.check` read 5 organs off the wrong PBPK state indices | `dS2dt` returns heart at 29–31, not in numeric order; kidney was checked against the heart ceiling, muscle against kidney, and so on. Invisible while all ceilings were 50.0, wrong as soon as one was set. |
| FES-1/FES-2 were fed raw `N_t` (~1e10) | Their `tumor_size` universe is `arange(-3, 12)`, i.e. log10(cells). scikit-fuzzy returns **zero** membership outside the universe, so tumour size contributed nothing and the controllers responded only to toxicity. Now fed `log10(N)`. |
| `pbpk[0..8]` unpacked unconditionally | Any cycle interval above 15 raised `IndexError`. |
| The repair loop had no iteration cap or dose floor | An unsatisfiable ceiling looped forever, each iteration a 35-ODE solve. |
| PBPK applied the dose twice on day 1 only when `int(dose) != 40` | `int()` truncates, so doses in [40, 41) were exempt: 39.99 gave ~1.8x the exposure of 40.00. Now applied uniformly. |
| Gompertz `ln(rho/N)` with N driven to zero | Short intervals eradicate the tumour, `N` crosses zero and NaN propagated through the solve and into the stored JSON. Growth is now evaluated at a positive floor. |
| Malformed UUID in `/compare?a=` | Returned a 500. Now ignored. |
| SciPy's LSODA is not reentrant | Runs and sweeps execute on worker threads, so a second simulation started while one was in flight raised `IntegratorConcurrencyError`. Integrations now take a shared lock (see `engine/_solver.py`). |
| `requirements.txt` was UTF-16LE and missing `httpx`, `whitenoise`, `gunicorn` | `pip install -r` failed, and the app could not boot from its own declared dependencies. |

### Performance

FES-1 and FES-2 rebuilt their entire scikit-fuzzy `ControlSystem` — 21 and 30 rules, all
membership functions, a networkx graph — on **every call**, from inside an ODE right-hand
side, and the `interval + 1` branch called both on every solver step only to discard the
result. Rule graphs are now built once at import.

A full run went from **~41 s to ~0.2 s**, and the result page from **2.6 MB to ~150 KB**
(charts are drawn client-side from JSON instead of ~40 inlined base64 matplotlib PNGs).

### Known model artifacts, not changed

Two PBPK equations look like transcription errors against the MATLAB but were left alone,
since changing model equations is a scientific decision:

- **EQ 7** (arterial RBC) uses `k_rbcplas` for the influx term where the symmetric EQ 2
  uses `k_plasrbc`.
- **EQ 15** (liver) reads `((F_li*C_art)+(F_g*C_gv)+(F_s*C_sv)-(C_liv*(F_g+F_s+F_li))/V_liv)`
  — the `/V_liv` binds only to the last term, so the inflow terms are missing their volume
  normalisation. Compare the correctly parenthesised EQ 3.

Also note the organ ceilings are **illustrative**. The source model defines none, and its
concentration units are not calibrated against clinical exposure limits.

## What the tool does beyond a single simulation

- **Regimen explorer** (`/explore`) sweeps cycle intervals against organ-limit profiles for
  one patient and plots tumour kill against peak toxicity. Regimens that nothing else beats
  on both axes are marked as the trade-off frontier. Nothing is ranked or recommended:
  which point is acceptable is a clinical judgement, not a model output.
- **Per-dose explanation** on every result: the state the controller saw, the linguistic
  terms it activated and with what membership, what FES-1 proposed, what FES-2 adjusted,
  whether the soft cap applied, and which organ ceiling cut the dose and by how much.
- **Patient-specific inputs.** Initial tumour burden is selectable in clinical terms
  (mass equivalents use the standard 1 g ~ 1e9 cells convention) and capped at 1e11, since
  the fuzzy universe ends at log10 = 11 and the Gompertz carrying capacity is 1e12.
  Body weight scales PBPK compartment volumes linearly and blood flows as BW^0.75, so a
  50 kg and a 120 kg patient no longer receive identical predicted organ exposures; at the
  70 kg reference weight every factor is exactly 1.0 and the published numbers are unchanged.
- **Run history and A/B comparison**, so a saved result can be reopened instantly.
