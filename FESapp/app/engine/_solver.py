"""Serialised access to SciPy's LSODA integrator.

The Fortran LSODA that ``solve_ivp(method="LSODA")`` wraps keeps global state and is not
reentrant: two threads integrating at the same time raise ``IntegratorConcurrencyError``.
Runs and regimen sweeps both execute on worker threads, so simultaneous simulations are
normal here - a user starting a second run, or a sweep, while one is still going.

The lock is taken per integration rather than per run. A single ODE solve is milliseconds
to tens of milliseconds, so a waiting request is never blocked for long, even behind a
regimen sweep that is running dozens of simulations.
"""
import threading

from scipy.integrate import solve_ivp

LSODA_LOCK = threading.Lock()


def solve(*args, **kwargs):
    """``scipy.integrate.solve_ivp`` guarded by the LSODA lock."""
    with LSODA_LOCK:
        return solve_ivp(*args, **kwargs)
