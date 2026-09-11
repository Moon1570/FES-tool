"""Runs simulations off the request thread and tracks their progress.

A full run is ~0.2 s, but an unlucky combination of organ limits can take several
seconds of repair iterations, and blocking the browser on a synchronous POST is how
the original hung. Runs execute on a worker thread; the browser polls for progress.

Progress lives in an in-process dict (fast, and the UI polls it several times a
second); the finished result is written to the database, which is what history and
comparison read. A single-process dev server is the deployment target, so an
in-memory registry is sufficient - a restart simply loses in-flight progress, and
completed runs survive in the database.
"""
import logging
import threading
import time
import traceback

from django.db import close_old_connections

from .engine import runner, serialize
from .models import Run

log = logging.getLogger(__name__)

_PROGRESS = {}
_LOCK = threading.Lock()


def set_progress(run_id, stage, fraction, message):
    with _LOCK:
        _PROGRESS[str(run_id)] = {
            "stage": stage,
            "progress": round(float(fraction), 4),
            "message": message,
        }


def get_progress(run_id):
    with _LOCK:
        return dict(_PROGRESS.get(str(run_id), {}))


def clear_progress(run_id):
    with _LOCK:
        _PROGRESS.pop(str(run_id), None)


def _execute(run_id):
    close_old_connections()
    started = time.time()
    try:
        run = Run.objects.get(pk=run_id)
        run.status = Run.Status.RUNNING
        run.save(update_fields=["status"])

        config = runner.RunConfig(
            patient_name=run.patient_name,
            weight_kg=run.weight_kg,
            interval_days=run.interval_days,
            horizon_days=run.horizon_days,
            organ_limits=run.organ_limits,
            n0=run.n0,
            log_tumor=run.log_tumor,
            allometric=run.allometric,
            legacy_dose_gate=run.legacy_dose_gate,
        )

        def progress(stage, fraction, message):
            set_progress(run_id, stage, fraction, message)

        set_progress(run_id, "start", 0.01, "Starting…")
        result = runner.run(config, progress=progress)

        run.result = serialize.serialize(result)
        run.duration_ms = int((time.time() - started) * 1000)
        run.status = Run.Status.DONE
        run.save(update_fields=["result", "duration_ms", "status"])
        set_progress(run_id, "done", 1.0, "Complete")
    except Exception as exc:  # surfaced to the user rather than swallowed
        log.exception("simulation failed for run %s", run_id)
        try:
            run = Run.objects.get(pk=run_id)
            run.status = Run.Status.ERROR
            run.error = f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}"
            run.duration_ms = int((time.time() - started) * 1000)
            run.save(update_fields=["status", "error", "duration_ms"])
        except Exception:
            log.exception("could not record failure for run %s", run_id)
        set_progress(run_id, "error", 1.0, str(exc))
    finally:
        close_old_connections()


def run_inline(run):
    """Execute ``run`` inside the current request and return when it is finished.

    For hosts that do not allow threads in web workers (PythonAnywhere runs uWSGI
    without --enable-threads, and users cannot change that). A run takes ~0.2 s.
    """
    _execute(str(run.id))
    clear_progress(run.id)


def start(run):
    """Kick off ``run`` on a worker thread and return immediately."""
    set_progress(run.id, "queued", 0.0, "Queued…")
    thread = threading.Thread(target=_execute, args=(str(run.id),),
                              name=f"fes-run-{run.id}", daemon=True)
    thread.start()
    return thread


def _execute_sweep(sweep_id):
    close_old_connections()
    started = time.time()
    from .engine import sweep as sweep_engine
    from .models import Sweep
    try:
        sw = Sweep.objects.get(pk=sweep_id)
        sw.status = "running"
        sw.save(update_fields=["status"])

        def progress(stage, fraction, message):
            set_progress(sweep_id, stage, fraction, message)

        set_progress(sweep_id, "start", 0.01, "Starting…")
        sw.result = sweep_engine.explore(
            weight_kg=sw.weight_kg, n0=sw.n0, intervals=sw.intervals,
            preset_keys=sw.preset_keys, horizon_days=sw.horizon_days,
            progress=progress,
        )
        sw.duration_ms = int((time.time() - started) * 1000)
        sw.status = "done"
        sw.save(update_fields=["result", "duration_ms", "status"])
        set_progress(sweep_id, "done", 1.0, "Complete")
    except Exception as exc:
        log.exception("sweep failed for %s", sweep_id)
        try:
            sw = Sweep.objects.get(pk=sweep_id)
            sw.status = "error"
            sw.error = f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}"
            sw.save(update_fields=["status", "error"])
        except Exception:
            log.exception("could not record sweep failure for %s", sweep_id)
        set_progress(sweep_id, "error", 1.0, str(exc))
    finally:
        close_old_connections()


def start_sweep(sweep):
    set_progress(sweep.id, "queued", 0.0, "Queued…")
    thread = threading.Thread(target=_execute_sweep, args=(str(sweep.id),),
                              name=f"fes-sweep-{sweep.id}", daemon=True)
    thread.start()
    return thread


def run_sweep_inline(sweep):
    """Inline counterpart of :func:`start_sweep`, for hosts without threads."""
    _execute_sweep(str(sweep.id))
    clear_progress(sweep.id)
