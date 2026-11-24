"""Shared paths and helpers of the Cry Atlas pipeline (cry_atlas/scripts/NN_*.py).  Run every script with `py -I` after sys.path.insert(0, scripts folder) (the scripts do this themselves)."""
import os, sys, json, time, threading, urllib.request, urllib.error, urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
CRY = os.path.dirname(HERE)                                  # cry_atlas
DATA = os.path.join(CRY, "data"); WORK = os.path.join(CRY, "work"); SITE = os.path.join(CRY, "site"); DOCS = os.path.join(CRY, "docs")
for d in (DATA, WORK, SITE, DOCS):
    os.makedirs(d, exist_ok=True)
# raw structure downloads (hundreds of MB, re-downloadable from the manifests): outside the project folder
RAW = os.environ.get("CRY_RAW") or os.path.join(os.environ.get("TEMP", os.path.expanduser("~")), "cry_raw")
os.makedirs(RAW, exist_ok=True)
UA = {"User-Agent": "atlas-cry-pipeline/1.0 (research use)"}
_lock = threading.Lock(); _last = {}


def get(url, binary=False, tries=5, timeout=300, gap=0.35, key="ncbi", data=None, headers=None):
    """GET / POST with a minimum gap between requests of one service (`key`), retries on transient errors; None on 400/404"""
    err = None
    for k in range(tries):
        with _lock:
            wait = gap - (time.time() - _last.get(key, 0.0))
            if wait > 0:
                time.sleep(wait)
            _last[key] = time.time()
        try:
            h = dict(UA); h.update(headers or {})
            req = urllib.request.Request(url, headers=h, data=data)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                b = r.read()
            return b if binary else b.decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None
            err = "HTTP %d" % e.code
            if e.code == 429:                                                  # rate limit: back off for real (the server asks clients to slow down)
                time.sleep(40 * (k + 1)); continue
        except Exception as e:
            err = str(e)
        time.sleep(3 * (k + 1))
    raise RuntimeError("%s: %s" % (url, err))


def read_fasta(path):
    out, h, s = [], None, []
    for ln in open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if h is not None:
                out.append((h, "".join(s)))
            h, s = ln[1:], []
        elif ln:
            s.append(ln.strip())
    if h is not None:
        out.append((h, "".join(s)))
    return out


def wsl(p):
    p = os.path.abspath(p)
    return "/mnt/" + p[0].lower() + p[2:].replace("\\", "/")


def jload(p, default=None):
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else default


def jdump(o, p, **kw):
    json.dump(o, open(p, "w", encoding="utf-8"), **kw)
