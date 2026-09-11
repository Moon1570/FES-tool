"""Views for the FES chemotherapy dose-scheduling tool.

All of the science lives in :mod:`app.engine`; these functions only marshal form
input, kick off background runs and render results.
"""
import csv
import json
import uuid

from django.conf import settings
from django.db.models import Q
from django.http import Http404, HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.views.decorators.http import require_GET, require_POST

from . import executor
from functools import lru_cache

from .engine import burden, runner, safety
from .engine import tumor as tumor_engine
from .engine import sweep as sweep_engine
from .models import Run, Sweep

MIN_WEIGHT, MAX_WEIGHT = 1.0, 400.0

# In public mode (FES_PUBLIC=1) the bounds are tighter. On a single-worker host one
# request blocks every other visitor, and the expensive inputs are tight organ limits
# (each dose is then cut 5% at a time, many times over) and large explorer grids.
# Measured worst cases under these bounds are about a second per run.
PUBLIC_MIN_INTERVAL = 7                 # weekly: the shortest common schedule
PUBLIC_MAX_HORIZON = 180
PUBLIC_LIMIT_FLOOR = 0.5                # fraction of the Standard preset per organ
PUBLIC_MAX_SWEEP_POINTS = 16
PUBLIC_SWEEP_INTERVALS = [7, 14, 21, 28]
REMEMBER_MAX = 50                       # run/sweep ids kept per visitor session


def _interval_bounds():
    return (PUBLIC_MIN_INTERVAL if settings.FES_PUBLIC else 1), 60


def _horizon_bounds():
    return 30, (PUBLIC_MAX_HORIZON if settings.FES_PUBLIC else 365)


def _limit_floor(organ_key):
    if settings.FES_PUBLIC:
        return PUBLIC_LIMIT_FLOOR * safety.STANDARD_LIMITS[organ_key]
    return 0.0001


def _dispatch_run(run):
    if settings.FES_RUN_MODE == "inline":
        executor.run_inline(run)
    else:
        executor.start(run)


def _dispatch_sweep(sweep):
    if settings.FES_RUN_MODE == "inline":
        executor.run_sweep_inline(sweep)
    else:
        executor.start_sweep(sweep)


# -- per-visitor visibility ---------------------------------------------------
# A public deployment has no accounts, so without this every visitor would see
# every other visitor's runs - including whatever they typed as the patient name.
# Locally (FES_PUBLIC unset) everything stays visible, as before.

def _remember(request, kind, obj_id):
    ids = [i for i in request.session.get(kind, []) if i != str(obj_id)]
    request.session[kind] = ([str(obj_id)] + ids)[:REMEMBER_MAX]


def _visible_runs(request):
    if not settings.FES_PUBLIC:
        return Run.objects.all()
    return Run.objects.filter(Q(is_demo=True) | Q(pk__in=request.session.get("runs", [])))


def _visible_sweeps(request):
    if not settings.FES_PUBLIC:
        return Sweep.objects.all()
    return Sweep.objects.filter(Q(is_demo=True) | Q(pk__in=request.session.get("sweeps", [])))


def _clamp(value, low, high):
    return max(low, min(high, value))


def _number(raw, default, low, high, cast=float):
    if raw is None or str(raw).strip() == "":
        return default
    try:
        return _clamp(cast(raw), low, high)
    except (TypeError, ValueError):
        return default


def _presets_for_template():
    return [
        {"key": key, "label": p["label"], "description": p["description"],
         "limits": p["limits"]}
        for key, p in safety.PRESETS.items()
    ]


def _burden_choices():
    """Burden presets annotated with mass and the fuzzy term each one activates."""
    out = []
    for p in burden.PRESETS:
        d = burden.describe(p["cells"])
        out.append(dict(p, mass=d["mass"], log10=d["log10"], fuzzy_term=d["fuzzy_term"]))
    return out


def new_run(request):
    return render(request, "run_setup.html", {
        "organs": safety.ORGANS,
        "presets": _presets_for_template(),
        "presets_json": json.dumps({k: p["limits"] for k, p in safety.PRESETS.items()}),
        "default_preset": safety.DEFAULT_PRESET,
        "burdens": _burden_choices(),
        "default_burden": burden.DEFAULT_PRESET,
        "defaults": {"weight": 70, "interval": 14, "horizon": 120},
        "bounds": {"interval": _interval_bounds(), "horizon": _horizon_bounds()},
        "public": settings.FES_PUBLIC,
        "recent": _visible_runs(request).filter(status=Run.Status.DONE)[:5],
    })


@require_POST
def create_run(request):
    """Validate the form, persist a Run and start it on a worker thread."""
    post = request.POST
    name = (post.get("name") or "").strip() or "John Doe"

    preset_key = post.get("preset") or safety.DEFAULT_PRESET
    preset = safety.PRESETS.get(preset_key, safety.PRESETS[safety.DEFAULT_PRESET])

    # Per-organ overrides fall back to the chosen preset rather than a bare 50.0.
    limits = {}
    for organ in safety.ORGANS:
        limits[organ.key] = _number(
            post.get(f"{organ.key}_max"), preset["limits"][organ.key],
            _limit_floor(organ.key), 10_000.0,
        )

    # Burden: a preset key, or "custom" with an explicit cell count.
    burden_key = post.get("burden") or burden.DEFAULT_PRESET
    if burden_key == "custom":
        n0 = _number(post.get("n0_custom"), float(burden.PRESET_BY_KEY[burden.DEFAULT_PRESET]["cells"]),
                     burden.MIN_CELLS, burden.MAX_CELLS)
    else:
        preset_b = burden.PRESET_BY_KEY.get(burden_key, burden.PRESET_BY_KEY[burden.DEFAULT_PRESET])
        n0 = preset_b["cells"]

    run = Run.objects.create(
        patient_name=name,
        n0=burden.clamp(n0),
        weight_kg=_number(post.get("weight"), 70.0, MIN_WEIGHT, MAX_WEIGHT),
        interval_days=_number(post.get("intervalTime"), 14, *_interval_bounds(), int),
        horizon_days=_number(post.get("horizon"), 120, *_horizon_bounds(), int),
        preset=preset_key if preset_key in safety.PRESETS else safety.DEFAULT_PRESET,
        organ_limits=limits,
        log_tumor=post.get("log_tumor", "1") == "1",
    )
    _remember(request, "runs", run.id)
    _dispatch_run(run)

    if request.headers.get("X-Requested-With") == "XMLHttpRequest":
        return JsonResponse({"id": str(run.id), "status": run.status})
    return redirect("result", run_id=run.id)


@require_GET
def run_status(request, run_id):
    run = get_object_or_404(Run, pk=run_id)
    payload = {"id": str(run.id), "status": run.status, "error": run.error}
    payload.update(executor.get_progress(run.id))
    return JsonResponse(payload)


@require_GET
def result(request, run_id):
    run = get_object_or_404(Run, pk=run_id)
    return render(request, "result.html", {
        "run": run,
        "organs": safety.ORGANS,
        # The same organ -> compartments map the safety check sums, so each chart
        # in the organ grid plots exactly the quantity compared with its limit.
        "organ_states_json": json.dumps({o.key: list(o.states) for o in safety.ORGANS}),
        "payload_json": json.dumps(run.result) if run.result else "null",
        "toxicity_limit": safety.LEGACY_DEFAULT_LIMIT,
    })


# Organs shown on the home page after the one(s) that limited the dose: the heart and
# kidneys are the organs clinicians watch most closely for chemotherapy damage.
CASE_STUDY_ORGANS = ["heart", "kidney", "liver", "brain"]


@lru_cache(maxsize=32)
def _untreated_log10(n0, interval_days, horizon_days, n_cycles):
    """The same tumour with every dose set to zero: the model's untreated course."""
    cfg = tumor_engine.TumorConfig(interval_days=interval_days, horizon_days=horizon_days, n0=n0)
    result = tumor_engine.simulate_replay([0.0] * n_cycles, cfg)
    return [round(float(v), 4) for v in result.log10N]


def _case_study(run):
    """Compact figures for the home page, all read from one stored demo run."""
    d = run.result
    by_key = {o["key"]: o for o in d["organs"]}
    states = {o.key: o.states for o in safety.ORGANS}

    picked = [o["key"] for o in d["organs"] if o["binding"]]
    for key in CASE_STUDY_ORGANS:
        if len(picked) >= 3:
            break
        if key not in picked:
            picked.append(key)

    organs = []
    for key in picked[:3]:
        o = by_key[key]
        totals = [[round(sum(c["states"][n][k] for n in states[key]), 4)
                   for k in range(len(d["pbpk"]["t_hours"]))]
                  for c in d["pbpk"]["cycles"]]
        organs.append({"key": key, "label": o["label"], "peak": o["peak"], "limit": o["limit"],
                       "binding": o["binding"], "series": totals})

    return {
        "run_id": str(run.id),
        "days": d["final"]["days"],
        "untreated": _untreated_log10(float(run.n0), run.interval_days, run.horizon_days,
                                      len(d["cycles"])),
        "treated": d["final"]["log10N"],
        "tox_planned": d["planned"]["toxicity"],
        "tox_final": d["final"]["toxicity"],
        "dose_planned": d["planned"]["dose_by_day"],
        "dose_final": d["final"]["dose_by_day"],
        "metrics": d["metrics"],
        "peak_tox_planned": max(v for v in d["planned"]["toxicity"] if v is not None),
        "hours": d["pbpk"]["t_hours"],
        "doses": [c["index"] for c in d["pbpk"]["cycles"]],
        "organs": organs,
    }


@require_GET
def features(request):
    """Static 'Key features' page. Its example links point at the seeded demo cases
    when they exist, and simply fall back to the setup page when they don't."""
    demo = Run.objects.filter(is_demo=True, status=Run.Status.DONE)
    example = demo.filter(patient_name="Case 01 — standard").first() or demo.first()
    # The card describes 14-day against 21-day cycles, so prefer that demo pair.
    other = (demo.filter(patient_name="Case 01 — 21-day cycles").first()
             or (demo.exclude(pk=example.pk).first() if example else None))
    return render(request, "features.html", {
        "example": example,
        "case": _case_study(example) if example and example.result else None,
        "compare_pair": (example, other) if example and other else None,
        "demo_sweep": Sweep.objects.filter(is_demo=True, status="done").first(),
    })


@require_GET
def history(request):
    return render(request, "history.html", {
        "runs": _visible_runs(request)[:100],
        "public": settings.FES_PUBLIC,
    })


def _lookup(runs, raw):
    """Resolve a run id from a query parameter, tolerating anything malformed.

    Django raises ValidationError (a 500) when a non-UUID string reaches a UUID
    primary-key filter, so a mangled or truncated link would crash this page.
    """
    if not raw:
        return None
    try:
        key = uuid.UUID(str(raw))
    except (ValueError, AttributeError, TypeError):
        return None
    return runs.filter(pk=key).first()


@require_GET
def compare(request):
    runs = _visible_runs(request).filter(status=Run.Status.DONE)
    left = _lookup(runs, request.GET.get("a"))
    right = _lookup(runs, request.GET.get("b"))
    return render(request, "compare.html", {
        "runs": runs[:100],
        "left": left,
        "right": right,
        "left_json": json.dumps(left.result) if left and left.result else "null",
        "right_json": json.dumps(right.result) if right and right.result else "null",
    })


@require_GET
def export_csv(request, run_id):
    """The dose schedule and daily curves, as a real download."""
    run = get_object_or_404(Run, pk=run_id)
    if not run.result:
        raise Http404("This run has no result yet.")
    data = run.result

    response = HttpResponse(content_type="text/csv")
    stem = "".join(ch if ch.isalnum() else "_" for ch in run.patient_name)[:40]
    response["Content-Disposition"] = f'attachment; filename="fes_{stem}_{run.id.hex[:8]}.csv"'
    writer = csv.writer(response)

    writer.writerow(["FES chemotherapy dose scheduling - research prototype, not for clinical use"])
    writer.writerow(["Patient", data["patient"]["name"]])
    writer.writerow(["Weight (kg)", data["patient"]["weight_kg"]])
    writer.writerow(["BSA (m2)", data["patient"]["bsa_m2"]])
    writer.writerow(["Cycle interval (days)", data["regimen"]["interval_days"]])
    writer.writerow(["Limit preset", run.preset_label])
    writer.writerow([])

    writer.writerow(["Cycle", "Day", "Planned dose", "Delivered dose",
                     "Reduction %", "Limiting organ", "Infeasible"])
    for c in data["cycles"]:
        writer.writerow([c["index"], c["day"], c["planned"], c["delivered"],
                         c["reduction_pct"], c["limiting_label"] or "", c["infeasible"]])
    writer.writerow([])

    writer.writerow(["Organ", "Limit", "Worst peak", "Utilisation", "Dose-limiting"])
    for o in data["organs"]:
        writer.writerow([o["label"], o["limit"], o["peak"], o["utilisation"], o["binding"]])
    writer.writerow([])

    final = data["final"]
    writer.writerow(["Day", "Dose", "Drug concentration", "log10(tumour cells)", "Toxicity"])
    for i, day in enumerate(final["days"]):
        writer.writerow([day, final["dose_by_day"][i], final["concentration"][i],
                         final["log10N"][i], final["toxicity"][i]])
    return response


# ---------------------------------------------------------------------------
# Regimen explorer
# ---------------------------------------------------------------------------

@require_GET
def explore(request):
    low, high = _interval_bounds()
    return render(request, "explore.html", {
        "intervals": [i for i in sweep_engine.DEFAULT_INTERVALS if low <= i <= high],
        "default_intervals": (PUBLIC_SWEEP_INTERVALS if settings.FES_PUBLIC
                              else sweep_engine.DEFAULT_INTERVALS),
        "max_points": PUBLIC_MAX_SWEEP_POINTS if settings.FES_PUBLIC else 0,
        "presets": _presets_for_template(),
        "default_presets": sweep_engine.DEFAULT_PRESETS,
        "burdens": _burden_choices(),
        "default_burden": burden.DEFAULT_PRESET,
        "recent": _visible_sweeps(request).filter(status="done")[:5],
    })


@require_POST
def create_sweep(request):
    post = request.POST
    intervals = [int(v) for v in post.getlist("intervals") if str(v).isdigit()]
    low, high = _interval_bounds()
    intervals = sorted({i for i in intervals if low <= i <= high})
    if not intervals:
        intervals = list(PUBLIC_SWEEP_INTERVALS if settings.FES_PUBLIC
                         else sweep_engine.DEFAULT_INTERVALS)

    preset_keys = [k for k in post.getlist("preset_keys") if k in safety.PRESETS]
    if not preset_keys:
        preset_keys = list(sweep_engine.DEFAULT_PRESETS)

    if settings.FES_PUBLIC:
        # Keep every chosen profile and trim intervals until the grid fits.
        per = max(1, PUBLIC_MAX_SWEEP_POINTS // len(preset_keys))
        intervals = intervals[:per]
        preset_keys = preset_keys[:PUBLIC_MAX_SWEEP_POINTS]

    burden_key = post.get("burden") or burden.DEFAULT_PRESET
    if burden_key == "custom":
        n0 = _number(post.get("n0_custom"), 1e10, burden.MIN_CELLS, burden.MAX_CELLS)
    else:
        n0 = burden.PRESET_BY_KEY.get(
            burden_key, burden.PRESET_BY_KEY[burden.DEFAULT_PRESET])["cells"]

    sweep = Sweep.objects.create(
        patient_name=(post.get("name") or "").strip() or "John Doe",
        weight_kg=_number(post.get("weight"), 70.0, MIN_WEIGHT, MAX_WEIGHT),
        n0=burden.clamp(n0),
        horizon_days=_number(post.get("horizon"), 120, *_horizon_bounds(), int),
        intervals=intervals, preset_keys=preset_keys,
    )
    _remember(request, "sweeps", sweep.id)
    _dispatch_sweep(sweep)
    return redirect("sweep_detail", sweep_id=sweep.id)


@require_GET
def sweep_detail(request, sweep_id):
    sweep = get_object_or_404(Sweep, pk=sweep_id)
    return render(request, "sweep.html", {
        "sweep": sweep,
        "payload_json": json.dumps(sweep.result) if sweep.result else "null",
        "burden": burden.describe(sweep.n0),
    })


@require_GET
def sweep_status(request, sweep_id):
    sweep = get_object_or_404(Sweep, pk=sweep_id)
    payload = {"id": str(sweep.id), "status": sweep.status, "error": sweep.error}
    payload.update(executor.get_progress(sweep.id))
    return JsonResponse(payload)


@require_POST
def run_from_sweep(request, sweep_id):
    """Turn one point of a sweep into a full, saved run."""
    sweep = get_object_or_404(Sweep, pk=sweep_id)
    interval = _number(request.POST.get("interval"), 14, *_interval_bounds(), int)
    preset_key = request.POST.get("preset")
    preset = safety.PRESETS.get(preset_key, safety.PRESETS[safety.DEFAULT_PRESET])

    run = Run.objects.create(
        patient_name=f"{sweep.patient_name} — {interval}-day, {preset['label']}",
        weight_kg=sweep.weight_kg, n0=sweep.n0,
        interval_days=interval, horizon_days=sweep.horizon_days,
        preset=preset_key if preset_key in safety.PRESETS else safety.DEFAULT_PRESET,
        organ_limits=preset["limits"],
    )
    _remember(request, "runs", run.id)
    _dispatch_run(run)
    return redirect("result", run_id=run.id)
