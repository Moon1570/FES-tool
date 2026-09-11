import uuid

from django.db import models

from .engine import burden, safety


class Run(models.Model):
    """One simulation run: its inputs, its status and its serialized result.

    Runs are persisted so a result can be reopened instantly (useful as a fallback
    during a live demo) and so two runs can be compared without re-simulating.
    """

    class Status(models.TextChoices):
        PENDING = "pending", "Pending"
        RUNNING = "running", "Running"
        DONE = "done", "Complete"
        ERROR = "error", "Failed"

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    created_at = models.DateTimeField(auto_now_add=True)

    patient_name = models.CharField(max_length=100, default="John Doe")
    weight_kg = models.FloatField(default=70.0)
    interval_days = models.PositiveIntegerField(default=14)
    horizon_days = models.PositiveIntegerField(default=120)
    n0 = models.FloatField(default=1e10, help_text="Initial tumour burden, in cells.")

    preset = models.CharField(max_length=32, default=safety.DEFAULT_PRESET)
    organ_limits = models.JSONField(default=dict)
    log_tumor = models.BooleanField(default=True)
    allometric = models.BooleanField(default=True)
    legacy_dose_gate = models.BooleanField(default=False)
    is_demo = models.BooleanField(
        default=False,
        help_text="Seeded example run. Visible to every visitor in public mode.")

    status = models.CharField(max_length=16, choices=Status.choices, default=Status.PENDING)
    error = models.TextField(blank=True, default="")
    result = models.JSONField(null=True, blank=True)
    duration_ms = models.PositiveIntegerField(null=True, blank=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self):
        return f"{self.patient_name} ({self.interval_days}-day) {self.created_at:%Y-%m-%d %H:%M}"

    @property
    def preset_label(self):
        return safety.PRESETS.get(self.preset, {}).get("label", self.preset)

    @property
    def burden_label(self):
        preset = burden.nearest_preset(self.n0)
        return preset["label"] if preset else f"{self.n0:.2e} cells"

    @property
    def summary(self):
        """Compact KPI dict for history and comparison listings."""
        if not self.result:
            return None
        m = self.result["metrics"]
        return {
            "log_reduction": m["log_reduction"],
            "peak_toxicity": m["peak_toxicity"],
            "cycles_adjusted": m["cycles_adjusted"],
            "cycles_total": m["cycles_total"],
            "total_delivered": m["total_delivered"],
        }


class Sweep(models.Model):
    """A grid of simulations across cycle intervals and organ-limit presets.

    A single run takes ~0.2 s, so exploring a few dozen regimens is a few seconds of
    work. The result is the efficacy/toxicity trade-off surface for one patient.
    """

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    created_at = models.DateTimeField(auto_now_add=True)

    patient_name = models.CharField(max_length=100, default="John Doe")
    weight_kg = models.FloatField(default=70.0)
    n0 = models.FloatField(default=1e10)
    horizon_days = models.PositiveIntegerField(default=120)

    intervals = models.JSONField(default=list)
    preset_keys = models.JSONField(default=list)
    is_demo = models.BooleanField(
        default=False,
        help_text="Seeded example sweep. Visible to every visitor in public mode.")

    status = models.CharField(max_length=16, choices=Run.Status.choices, default=Run.Status.PENDING)
    error = models.TextField(blank=True, default="")
    result = models.JSONField(null=True, blank=True)
    duration_ms = models.PositiveIntegerField(null=True, blank=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self):
        return f"Sweep for {self.patient_name} ({len(self.intervals)}x{len(self.preset_keys)})"

    @property
    def point_count(self):
        return len(self.intervals) * len(self.preset_keys)
