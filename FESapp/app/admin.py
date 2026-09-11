from django.contrib import admin

from .models import Run


@admin.register(Run)
class RunAdmin(admin.ModelAdmin):
    list_display = ("patient_name", "weight_kg", "interval_days", "preset", "status", "created_at")
    list_filter = ("status", "preset")
    readonly_fields = ("id", "created_at", "result", "duration_ms")
