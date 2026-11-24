"""Fetch structure files of the Cry rows (only alphafold.ebi.ac.uk/files and files.rcsb.org are contacted; files are untrusted data and are never executed; read with py -I).
 AFDB : work/afdb/map.json -> RAW/afdb/<entry>-model_v<N>.pdb   (exact-sequence models only)
 RCSB : work/pdb/matches.json (entries with an entity >= 90 % identical to a row) -> RAW/pdb/<ID>.cif
Manifest (url, bytes, sha256) -> work/structures_manifest.tsv ; resumable.  Run: py -I 12_fetch_structures.py [afdb|pdb|both]"""
import os, sys, re, json, csv, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
what = sys.argv[1] if len(sys.argv) > 1 else "both"
os.makedirs(os.path.join(RAW, "afdb"), exist_ok=True); os.makedirs(os.path.join(RAW, "pdb"), exist_ok=True)
MAN = os.path.join(WORK, "structures_manifest.tsv")
man = {}
if os.path.exists(MAN):
    for r in csv.DictReader(open(MAN, encoding="utf-8"), delimiter="\t"): man[r["file"]] = r


def save_manifest():
    with open(MAN, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["file", "source", "url", "bytes", "sha256"])
        for k, r in sorted(man.items()): w.writerow([k, r["source"], r["url"], r["bytes"], r["sha256"]])


def fetch(url, dest, source, key, gap):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        rel = os.path.relpath(dest, RAW).replace("\\", "/")
        if rel not in man:                                                  # downloaded by an interrupted earlier run: record it
            b = open(dest, "rb").read(); man[rel] = {"file": rel, "source": source, "url": url, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
        return True
    b = get(url, binary=True, key=key, gap=gap)
    if not b or (b.find(b"_atom_site.") < 0 and b.find(b"\nATOM  ") < 0):                      # neither an mmCIF atom_site loop nor PDB ATOM records
        return False
    open(dest, "wb").write(b)
    man[os.path.relpath(dest, RAW).replace("\\", "/")] = {"file": os.path.relpath(dest, RAW).replace("\\", "/"), "source": source, "url": url, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
    return True


if what in ("afdb", "both"):
    mp = jload(os.path.join(WORK, "afdb", "map.json")); ok = bad = 0
    ents = sorted({v["afdb"]["entry"] + "|" + v["afdb"]["pdbUrl"] for v in mp.values() if v.get("afdb")})
    print("AFDB models to fetch:", len(ents), flush=True)
    from concurrent.futures import ThreadPoolExecutor
    def one(e):
        entry, url = e.split("|"); dest = os.path.join(RAW, "afdb", os.path.basename(url))
        try: return url, fetch(url, dest, "AlphaFold DB", "afdb", 0.5)
        except Exception as ex: return url, False
    with ThreadPoolExecutor(2) as ex:
        for k, (url, good) in enumerate(ex.map(one, ents)):
            if good: ok += 1
            else: bad += 1; print("  failed", url, flush=True)
            if (k + 1) % 50 == 0: print("  AFDB %d of %d" % (k + 1, len(ents)), flush=True)
    save_manifest()
    save_manifest(); print("AFDB done: ok", ok, "failed", bad, flush=True)

if what in ("pdb", "both"):
    ids = set(jload(os.path.join(WORK, "pdb", "matches.json"), {}))                 # entries with an entity >= 90 % identical to a row (08c_pdb_entities.py)
    print("PDB entries to fetch:", len(ids), flush=True); ok = bad = 0
    for k, i in enumerate(sorted(ids)):
        dest = os.path.join(RAW, "pdb", i + ".cif")
        if fetch("https://files.rcsb.org/download/%s.cif" % i, dest, "RCSB PDB", "rcsb", 0.3): ok += 1
        else: bad += 1; print("  failed", i, flush=True)
    save_manifest(); print("PDB done: ok", ok, "failed", bad, flush=True)
