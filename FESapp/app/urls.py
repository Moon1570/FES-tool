from django.urls import path

from . import views

urlpatterns = [
    path("", views.home, name="home"),
    path("runs/new", views.create_run, name="create_run"),
    path("runs/<uuid:run_id>", views.result, name="result"),
    path("runs/<uuid:run_id>/status", views.run_status, name="run_status"),
    path("runs/<uuid:run_id>/export.csv", views.export_csv, name="export_csv"),
    path("explore", views.explore, name="explore"),
    path("explore/new", views.create_sweep, name="create_sweep"),
    path("explore/<uuid:sweep_id>", views.sweep_detail, name="sweep_detail"),
    path("explore/<uuid:sweep_id>/status", views.sweep_status, name="sweep_status"),
    path("explore/<uuid:sweep_id>/run", views.run_from_sweep, name="run_from_sweep"),
    path("history", views.history, name="history"),
    path("compare", views.compare, name="compare"),
]
