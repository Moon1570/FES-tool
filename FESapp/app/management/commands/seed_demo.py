"""Populate the database with a small set of presentable runs.

Useful two ways: History and Compare have something in them straight away, and a
saved run is a fallback if a live run misbehaves during a presentation.

    python manage.py seed_demo --reset
"""
import time

from django.core.management.base import BaseCommand

from app.engine import runner, safety, serialize, sweep as sweep_engine
from app.models import Run, Sweep

SCENARIOS = [
    ("Case 01 — standard", 70.0, 14, "standard"),
    ("Case 01 — 21-day cycles", 70.0, 21, "standard"),
    ("Case 02 — reduced renal", 82.0, 14, "renal"),
]


class Command(BaseCommand):
    help = "Create a set of demo runs."

    def add_arguments(self, parser):
        parser.add_argument("--reset", action="store_true",
                            help="Delete every existing run first.")

    def handle(self, *args, **options):
        if options["reset"]:
            runs, _ = Run.objects.all().delete()
            sweeps, _ = Sweep.objects.all().delete()
            self.stdout.write(f"Deleted {runs} run row(s) and {sweeps} sweep row(s).")
        else:
            # Without --reset, replace the demo cases instead of duplicating them,
            # and leave visitors' own runs alone.
            Run.objects.filter(is_demo=True).delete()
            Sweep.objects.filter(is_demo=True).delete()

        for name, weight, interval, preset in SCENARIOS:
            limits = safety.PRESETS[preset]["limits"]
            started = time.time()
            result = runner.run(runner.RunConfig(
                patient_name=name, weight_kg=weight, interval_days=interval,
                organ_limits=limits,
            ))
            payload = serialize.serialize(result)
            Run.objects.create(
                patient_name=name, weight_kg=weight, interval_days=interval,
                horizon_days=120, preset=preset, organ_limits=limits,
                status=Run.Status.DONE, result=payload, is_demo=True,
                duration_ms=int((time.time() - started) * 1000),
            )
            m = payload["metrics"]
            self.stdout.write(
                f"  {name:28} log-kill {m['log_reduction']:5.2f}  "
                f"peak tox {m['peak_toxicity']:5.1f}  "
                f"{m['cycles_adjusted']}/{m['cycles_total']} adjusted"
            )
        intervals = [7, 10, 14, 17, 21, 24, 28]
        presets = ["legacy", "standard", "renal", "cardiac"]
        started = time.time()
        payload = sweep_engine.explore(
            weight_kg=70.0, n0=1e10, intervals=intervals,
            preset_keys=presets, horizon_days=120,
        )
        Sweep.objects.create(
            patient_name="Case 01 — regimen sweep", weight_kg=70.0, n0=1e10,
            horizon_days=120, intervals=intervals, preset_keys=presets,
            status=Run.Status.DONE, result=payload, is_demo=True,
            duration_ms=int((time.time() - started) * 1000),
        )
        front = sum(1 for p in payload["points"] if p["on_front"])
        self.stdout.write(
            f"  {'Case 01 — regimen sweep':28} {len(payload['points'])} regimens, "
            f"{front} on the trade-off frontier"
        )
        self.stdout.write(self.style.SUCCESS(
            f"Seeded {len(SCENARIOS)} runs and 1 regimen sweep."))
