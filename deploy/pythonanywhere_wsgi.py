# PythonAnywhere WSGI configuration for FES-tool.
#
# Paste this whole file over the contents of the WSGI configuration file linked from
# your Web tab (/var/www/<username>_pythonanywhere_com_wsgi.py), then change the two
# values marked CHANGE-ME. See DEPLOY.md, step 6.
import os
import sys

USERNAME = "CHANGE-ME"      # your PythonAnywhere username, e.g. "fesdemo"
SECRET_KEY = "CHANGE-ME"    # the key you generated in DEPLOY.md step 4

# Accounts created on eu.pythonanywhere.com live at <username>.eu.pythonanywhere.com.
HOST = f"{USERNAME}.pythonanywhere.com"

if "CHANGE-ME" in (USERNAME, SECRET_KEY):
    raise RuntimeError("Edit USERNAME and SECRET_KEY at the top of the WSGI file (DEPLOY.md step 6).")

project = f"/home/{USERNAME}/FES-tool/FESapp"
if project not in sys.path:
    sys.path.insert(0, project)

os.environ["DJANGO_SETTINGS_MODULE"] = "FESapp.settings"
os.environ["DJANGO_DEBUG"] = "0"
os.environ["DJANGO_ALLOWED_HOSTS"] = HOST
os.environ["DJANGO_SECRET_KEY"] = SECRET_KEY
os.environ["FES_RUN_MODE"] = "inline"   # PythonAnywhere does not allow threads in web apps
os.environ["FES_PUBLIC"] = "1"          # per-visitor history and tighter input limits

from django.core.wsgi import get_wsgi_application  # noqa: E402

application = get_wsgi_application()
