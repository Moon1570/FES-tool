"""Regression check: new engine vs the recorded output of the original views.calc.

Legacy mode (raw tumour input, 50.0 organ ceilings) must reproduce
``tests/baseline_original.json`` exactly. Then the log10 fix is enabled and the
difference it makes is reported, so the port and the science change stay separable.
"""
import json
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from app.engine import runner, safety  # noqa: E402

BASELINE = os.path.join(os.path.dirname(HERE), "tests", "baseline_original.json")
TOL = 1e-9


def run_case(weight, interval, log_tumor, limits=None, legacy_gate=False,
             allometric=None):
    # Legacy checks must reproduce the published 70 kg reference physiology exactly,
    # so allometric scaling defaults off whenever the legacy dose gate is on.
    if allometric is None:
        allometric = not legacy_gate
    cfg = runner.RunConfig(
        patient_name="regress", weight_kg=weight, interval_days=interval,
        organ_limits=limits or {}, log_tumor=log_tumor,
        legacy_dose_gate=legacy_gate, allometric=allometric,
    )
    t0 = time.time()
    res = runner.run(cfg)
    return res, time.time() - t0


def main():
    with open(BASELINE) as fh:
        base = json.load(fh)

    failures = []
    print("=" * 78)
    print("1. LEGACY MODE must reproduce the original's FIRST (uncontaminated) run")
    print("=" * 78)
    expected = base["default_70kg_i14"]
    res, secs = run_case(70.0, 14, log_tumor=False, legacy_gate=True)
    got = [c.planned_dose for c in res.cycles]
    want = expected["doses"]
    ok_doses = len(got) == len(want) and all(abs(a - b) < TOL for a, b in zip(got, want))
    ok_bsa = abs(res.bsa_scaled - float(expected["BSA"])) < 1e-9
    ok_n84 = abs(res.planned.N[84] - float(expected["cell84"])) < 1e-9
    if not (ok_doses and ok_bsa and ok_n84):
        failures.append("legacy-equivalence")
    print(f"  [{'PASS' if ok_doses and ok_bsa and ok_n84 else 'FAIL'}] "
          f"default_70kg_i14  {secs:5.2f}s  "
          f"doses={'ok' if ok_doses else 'MISMATCH'} "
          f"BSA={'ok' if ok_bsa else 'MISMATCH'} "
          f"N[84]={'ok' if ok_n84 else 'MISMATCH'}")
    if not ok_doses:
        print(f"          got  {[round(x, 6) for x in got]}")
        print(f"          want {[round(x, 6) for x in want]}")
    speedup = expected["elapsed_sec"] / max(secs, 1e-6)
    print(f"  original: {expected['elapsed_sec']}s and "
          f"{expected['response_bytes'] / 1e6:.1f} MB per page  ->  {speedup:.0f}x faster")

    print()
    print("=" * 78)
    print("1b. The recorded 2nd run PROVES the module-global bug")
    print("=" * 78)
    # The original ran 90 kg immediately after 70 kg in one process. Its recorded
    # doses are byte-identical to the 70 kg run because the global `dose` list was
    # memoised and never reset, so weight had no effect at all. A correct engine
    # therefore must NOT reproduce it.
    leaked = base["second_run_90kg_i14"]["doses"] == base["default_70kg_i14"]["doses"]
    res90, _ = run_case(90.0, 14, log_tumor=False, legacy_gate=True)
    got90 = [c.planned_dose for c in res90.cycles]
    differs = any(abs(a - b) > 1e-6 for a, b in zip(got90, want))
    if leaked and differs:
        print("  [PASS] original leaked the 70 kg schedule into the 90 kg run;")
        print(f"         new engine gives 90 kg -> {[round(x, 2) for x in got90]}")
    else:
        failures.append("global-state-evidence")
        print(f"  [FAIL] leaked={leaked} differs={differs}")

    print()
    print("=" * 78)
    print("2. Two different patients must NOT share a schedule (the global-state bug)")
    print("=" * 78)
    a, _ = run_case(70.0, 14, log_tumor=True)
    b, _ = run_case(90.0, 14, log_tumor=True)
    da = [round(c.planned_dose, 4) for c in a.cycles]
    db = [round(c.planned_dose, 4) for c in b.cycles]
    if da == db:
        failures.append("patient-independence")
        print(f"  [FAIL] identical schedules: {da}")
    else:
        print(f"  [PASS] 70 kg -> {da}")
        print(f"         90 kg -> {db}")

    print()
    print("=" * 78)
    print("3. Dose interval must not crash (original raised IndexError for >15)")
    print("=" * 78)
    for interval in (7, 14, 21, 28, 40):
        try:
            res, secs = run_case(70.0, interval, log_tumor=True)
            print(f"  [PASS] interval={interval:<3} {len(res.cycles)} cycles, {secs:5.2f}s")
        except Exception as exc:
            failures.append(f"interval-{interval}")
            print(f"  [FAIL] interval={interval:<3} {type(exc).__name__}: {exc}")

    print()
    print("=" * 78)
    print("4. Unsatisfiable organ limit must terminate, not hang")
    print("=" * 78)
    t0 = time.time()
    res, _ = run_case(70.0, 14, log_tumor=True, limits={"kidney": 0.0001})
    secs = time.time() - t0
    infeasible = [c.index for c in res.cycles if c.infeasible]
    if infeasible and secs < 120:
        print(f"  [PASS] terminated in {secs:.1f}s; cycles flagged infeasible: {infeasible}")
    else:
        failures.append("infeasible-guard")
        print(f"  [FAIL] {secs:.1f}s, infeasible={infeasible}")

    print()
    print("=" * 78)
    print("5. Effect of the log10 fix (science change, reported not asserted)")
    print("=" * 78)
    old, _ = run_case(70.0, 14, log_tumor=False)
    new, _ = run_case(70.0, 14, log_tumor=True)
    print(f"  raw N_t  : {[round(c.planned_dose, 2) for c in old.cycles]}")
    print(f"  log10(N) : {[round(c.planned_dose, 2) for c in new.cycles]}")
    print(f"  log kill : raw {old.metrics['log_reduction']:.2f} -> "
          f"log10 {new.metrics['log_reduction']:.2f}")
    print(f"  peak tox : raw {old.metrics['peak_toxicity']:.2f} -> "
          f"log10 {new.metrics['peak_toxicity']:.2f}")

    print()
    print("=" * 78)
    print("5b. Body weight must change PBPK exposure (allometric scaling)")
    print("=" * 78)
    from app.engine import pbpk as _pbpk  # noqa: E402
    ref = _pbpk.simulate_pbpk(40.0)
    at70 = _pbpk.simulate_pbpk(40.0, phys=_pbpk.Physiology(70.0))
    identical = float(abs(ref.y - at70.y).max()) == 0.0
    kidney = safety.ORGAN_BY_KEY["kidney"]
    peaks = [safety.organ_peak(_pbpk.simulate_pbpk(40.0, phys=_pbpk.Physiology(w)).y, kidney)
             for w in (50, 70, 90, 120)]
    decreasing = all(b < a for a, b in zip(peaks, peaks[1:]))
    if identical and decreasing:
        print("  [PASS] 70 kg identical to the published reference (diff 0.0)")
        print("         kidney peak 50/70/90/120 kg: " +
              "  ".join(f"{p:.2f}" for p in peaks) + "  (falls with weight, as expected)")
    else:
        failures.append("allometric-scaling")
        print(f"  [FAIL] identical_at_70={identical} decreasing={decreasing} {peaks}")

    print()
    print("=" * 78)
    print("5c. Concurrent simulations must not collide (LSODA is not reentrant)")
    print("=" * 78)
    import threading  # noqa: E402
    errors, got = [], []

    def _go(w):
        try:
            got.append((w, round(run_case(w, 14, True, safety.STANDARD_LIMITS)[0]
                                 .metrics["log_reduction"], 3)))
        except Exception as exc:
            errors.append((w, type(exc).__name__, str(exc)[:80]))

    threads = [threading.Thread(target=_go, args=(w,)) for w in (55, 65, 70, 80, 90, 100)]
    t0 = time.time()
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    if errors or len(got) != len(threads):
        failures.append("solver-concurrency")
        print(f"  [FAIL] {errors}")
    else:
        print(f"  [PASS] {len(threads)} concurrent runs, no collisions, {time.time() - t0:.2f}s")
        print("         " + "  ".join(f"{w}kg:{v}" for w, v in sorted(got)))

    print()
    print("=" * 78)
    print("6. Dose -> exposure must be monotone (the int(dose)!=40 discontinuity)")
    print("=" * 78)
    from app.engine import pbpk  # noqa: E402
    kidney = safety.ORGAN_BY_KEY["kidney"]
    probes = [39.99, 40.0, 40.5, 40.99, 41.0, 42.5]
    peaks = [safety.organ_peak(pbpk.simulate_pbpk(d).y, kidney) for d in probes]
    monotone = all(b >= a for a, b in zip(peaks, peaks[1:]))
    jump = max(abs(b - a) / a for a, b in zip(peaks, peaks[1:]))
    if monotone and jump < 0.10:
        print(f"  [PASS] monotone, largest step {jump * 100:.1f}%")
    else:
        failures.append("dose-monotonicity")
        print(f"  [FAIL] monotone={monotone} largest step {jump * 100:.1f}%")
    print("         " + "  ".join(f"{d}->{p:.2f}" for d, p in zip(probes, peaks)))

    print()
    print("=" * 78)
    print("7. Organ-limit presets")
    print("=" * 78)
    for name, preset in safety.PRESETS.items():
        res, _ = run_case(70.0, 14, log_tumor=True, limits=preset["limits"])
        m = res.metrics
        organs = sorted({c.limiting_label for c in res.cycles if c.adjusted} - {None})
        print(f"  {name:9} {m['cycles_adjusted']}/{m['cycles_total']} cycles adjusted  "
              f"log-kill {m['log_reduction']:5.2f}  peak tox {m['peak_toxicity']:5.1f}  "
              f"limited by: {', '.join(organs) or '-'}")

    print()
    if failures:
        print(f"FAILURES: {failures}")
        return 1
    print("ALL REGRESSION CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
