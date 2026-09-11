# Deploying FES-tool for the fair

Free, always on, no cold start. About **45 minutes** end to end.

You will end up with three addresses:

| What | Address | Purpose |
|---|---|---|
| Live app | `https://YOURNAME.pythonanywhere.com` | The tool itself, on PythonAnywhere's free tier |
| Landing page | `https://moon1570.github.io/FES-tool/` | **What the QR code points to.** Two buttons: live demo and video |
| Video | `https://youtu.be/...` | Optional two-minute walkthrough |

Wherever you see **`YOURNAME`** below, use the PythonAnywhere username you pick in step 1.

> **Why the QR points at the landing page, not the app:** the printed poster can't be
> changed, but the landing page can. If the app moves host, lapses, or you add the video
> later, you edit one link and the QR keeps working.

---

## Part 0 — Push the code (on your laptop, 2 min)

Review what will be committed, then push:

```bash
cd ~/Research/chemo_drug_scheduling/FES-tool
git status                  # the deleted assests/ files are intended
git add -A
git commit -m "Revamp app and prepare PythonAnywhere deployment"
git push origin master
```

`.venv/`, `db.sqlite3` and `staticfiles/` are git-ignored, so they will not be uploaded.

---

## Part A — The live app on PythonAnywhere (~30 min)

### 1. Create a free account

Go to <https://www.pythonanywhere.com> → **Pricing & signup** → create the **free** account (currently called *Beginner*).

- The **username becomes your web address**, so keep it short and lowercase
  (for example `fesdemo` → `fesdemo.pythonanywhere.com`).
- Sign up on **www**.pythonanywhere.com, not eu.pythonanywhere.com. EU accounts get a
  different address (`YOURNAME.eu.pythonanywhere.com`); if you already have one, change
  the `HOST` line in the WSGI file in step 6.
- Confirm your email address.

### 2. Open a Bash console

**Dashboard** → **Consoles** → **Bash**. The remaining commands in Part A run here.

### 3. Download the code and install

```bash
git clone --depth 1 https://github.com/Moon1570/FES-tool.git
python3.11 -m venv ~/.virtualenvs/fes
source ~/.virtualenvs/fes/bin/activate
python -c "import _posixsubprocess, sys; print('venv ok', sys.version.split()[0])"

# numpy and scipy first, then drop the test suites they bundle, then the rest
pip install --no-cache-dir numpy==1.23.5 scipy==1.9.3
find ~/.virtualenvs/fes/lib/python3.11/site-packages -type d -name tests -prune -exec rm -rf {} +
pip install --no-cache-dir -r ~/FES-tool/FESapp/requirements.txt
find ~/.virtualenvs/fes/lib/python3.11/site-packages -type d -name tests -prune -exec rm -rf {} +

# must print: engine ok, 9 cycles, 12.27 log kill
cd ~/FES-tool/FESapp
python -c "from app.engine import runner; r = runner.run(runner.RunConfig()); print('engine ok,', len(r.cycles), 'cycles,', round(r.metrics['log_reduction'],2), 'log kill')"
```

Create the virtualenv with Python's own `venv` module, from the `python3.11` on your
PATH, so the interpreter and its standard library come from the same installation.
**Don't use `mkvirtualenv --python=/usr/bin/python3.11`**: on current PythonAnywhere
images that mixes two installations, and pip then fails with
`No module named '_posixsubprocess'`. The `venv ok` line must print before you install.
If `-m venv` reports that *ensurepip is not available*, use
`virtualenv -p "$(which python3.11)" ~/.virtualenvs/fes` instead.

Takes a few minutes. The finished install is about **380 MB of the 512 MB** free quota, but pip briefly needs more than that while it unpacks. Installing it in this order, and deleting the packages' bundled test suites (about 50 MB the app never uses), keeps the peak under the limit.

- Keep `--no-cache-dir`. Without it pip keeps a second copy of every download and you
  run out of space.
- Use **Python 3.11** (3.10 also works). **Not 3.12 or later** — the pinned numpy and
  scipy versions have no builds for it.
- The prompt should now start with `(fes)`. If you open a new console later, run
  `workon fes` first.

### 4. Database, demo cases, static files, secret key

```bash
cd ~/FES-tool/FESapp
python manage.py migrate
python manage.py seed_demo --reset
python manage.py collectstatic --noinput
python -c "import secrets; print(secrets.token_urlsafe(50))"
```

`seed_demo` creates the three demo cases and the 28-schedule comparison that every visitor
sees. **Copy the long random string** the last command prints; you need it in step 6.

### 5. Create the web app

**Web** tab → **Add a new web app** → Next → choose **Manual configuration**
(*not* "Django" — that option creates a new, empty project) → **Python 3.11** → Next.

### 6. Configure it

Still on the **Web** tab, fill in these sections:

**Code**
- Source code: `/home/YOURNAME/FES-tool/FESapp`
- Working directory: `/home/YOURNAME/FES-tool/FESapp`
- WSGI configuration file: click the link. Delete everything in that file and paste in
  the whole of [`deploy/pythonanywhere_wsgi.py`](deploy/pythonanywhere_wsgi.py). Then edit
  the two lines at the top:
  ```python
  USERNAME = "YOURNAME"
  SECRET_KEY = "the-string-from-step-4"
  ```
  Click **Save**.

**Virtualenv**
- `/home/YOURNAME/.virtualenvs/fes`

**Static files** — add one row:
| URL | Directory |
|---|---|
| `/static/` | `/home/YOURNAME/FES-tool/FESapp/staticfiles` |

**Security**
- **Force HTTPS**: turn it **on**. Forms are rejected over plain http.

Click the green **Reload YOURNAME.pythonanywhere.com** button at the top.

### 7. Test it

Open `https://YOURNAME.pythonanywhere.com` on your laptop **and on your phone**:

- [ ] The page is styled (blue header, cards). If it is plain text, see *Troubleshooting*.
- [ ] **Run simulation** with the defaults opens a finished result straight away.
- [ ] **History** lists the three "Case 0x" demo cases.
- [ ] **Explore** lists "Case 01 — regimen sweep", and it opens.
- [ ] Open the site in a private window: your test run should **not** appear in its History.

### 8. Stop it from expiring

Free web apps switch off after **one month** unless extended.
**Web** tab → the **Run until … from today** button (near the top). Put a monthly reminder in your calendar.

---

## Part B — The landing page on GitHub Pages (~5 min)

### 9. Put your address into the landing page

On your laptop, open `docs/index.html` and replace:

```
https://REPLACE-WITH-YOUR-USERNAME.pythonanywhere.com/
```

with your real address. Leave the YouTube link as it is for now — **the video button stays
hidden until its link is filled in**, so the page never shows a broken button.

```bash
git add docs/index.html
git commit -m "Point landing page at the live app"
git push origin master
```

### 10. Turn on GitHub Pages

<https://github.com/Moon1570/FES-tool> → **Settings** → **Pages** →
*Build and deployment*: Source **Deploy from a branch** → Branch **master**, folder
**/docs** → **Save**.

Within a minute or two it is live at **<https://moon1570.github.io/FES-tool/>**.

### 11. Test it on your phone

Open that address and tap **Try the live demo**.

---

## Part C — The video (optional, whenever you're ready)

12. Record a two-minute screen capture. A running order that works:
    *New run → Run simulation → treatment schedule (four doses lowered, lung is
    dose-limiting) → organ safety panel → Compare 14-day against 21-day cycles.*
13. Upload to YouTube as **Unlisted** or **Public** (not Private — nobody else can open a
    Private video). Copy the share link.
14. In `docs/index.html`, replace `https://youtu.be/REPLACE-WITH-VIDEO-ID` with your link,
    then commit and push. The **Watch the demo video** button appears within a minute.
    Nothing on the poster needs to change.

---

## Part D — The QR code

It should encode **`https://moon1570.github.io/FES-tool/`**. It is generated as part of the
poster rebuild.

If you make one yourself, use a **static** QR code. Many free "dynamic QR" sites encode a
redirect through their own servers and stop working when a trial ends — after the poster is
printed.

---

## The day before the fair

- [ ] Web tab → **Run until … from today**, so it can't lapse during the fair.
- [ ] Scan the printed QR with two different phones and run one simulation each.
- [ ] Laptop fallback works: `cd FESapp && ../.venv/bin/python manage.py runserver`.
- [ ] The video link is in the landing page, if you made one.

---

## What's different on the public site

It runs with `FES_PUBLIC=1` (set in the WSGI file). Your laptop is unaffected.

| | Laptop | Public site |
|---|---|---|
| History / Compare / Explore | everything | **each visitor's own runs, plus the demo cases** |
| Cycle interval | 1–60 days | 7–60 days |
| Treatment horizon | up to 365 days | up to 180 days |
| Organ limits | any value | no lower than half the Standard preset |
| Explorer | any grid | up to 16 schedules (the seeded 28-schedule demo is still shown) |

The free tier handles **one request at a time**. A normal run takes a fraction of a second,
so this is fine for fair traffic. The Explorer takes a few seconds, and anyone else
loading a page at that moment waits for it.

The limits exist because a single visitor with very tight organ limits could otherwise
keep the server busy for over ten seconds.

---

## Updating the live site later

In a PythonAnywhere Bash console:

```bash
cd ~/FES-tool && git pull
workon fes
cd FESapp
python manage.py migrate
python manage.py collectstatic --noinput
```

Then **Web** tab → **Reload**.

---

## Troubleshooting

Start with **Web** tab → **Error log** (the newest lines are at the bottom).

| Symptom | Cause and fix |
|---|---|
| "Something went wrong :-(" | Read the error log. If it says *Edit USERNAME and SECRET_KEY*, you haven't filled those in (step 6). |
| **400 Bad Request** | `USERNAME` in the WSGI file doesn't match your address — or it's an EU account, see step 1. |
| **403 Forbidden** when you press Run | You're on `http://`. Turn on **Force HTTPS** (step 6) and use `https://`. |
| Page loads with no styling | The static files row is wrong, or `collectstatic` wasn't run (steps 4 and 6). |
| `No module named 'django'` in the error log | The virtualenv path is wrong (step 6), or the install failed (step 3). |
| `No module named '_posixsubprocess'` when running pip | The virtualenv's interpreter and standard library come from different Python installations. `deactivate`, `rm -rf ~/.virtualenvs/fes`, and recreate it as in step 3. |
| **Disk quota exceeded** during `pip install` | Run the `find … -name tests …` line from step 3 and `rm -rf ~/.cache/pip`, then run the same `pip install --no-cache-dir -r …` again; it only installs what is missing. Then run `pip install --no-cache-dir --force-reinstall --no-deps matplotlib==3.6.2` and `pip check`, in case the last package was cut short. |
| The console gets very slow | Free accounts get 100 CPU-seconds a day for consoles. The allowance resets daily and **does not affect the website**. |
| The site stopped working after a few weeks | It expired. Click **Run until … from today** (step 8). |
