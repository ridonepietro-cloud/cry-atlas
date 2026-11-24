"""Quick test of cry_pack.build on a handful of rows whose RING edge files exist (first-pass frame work/aln/seq.fa); writes site/data/cry_test/ (not used by the page).  Run: py -I 19t_test_pack.py [n]"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
import cry_pack as CP
n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
rows = jload(os.path.join(WORK, "rows.json")); first = dict(CC.read_fasta(os.path.join(WORK, "aln", "seq.fa"))); mods = jload(os.path.join(WORK, "modules_seq.json"))["modules"]
sel = []
for r in rows:
    ms = [m for m in r["models"] if m.get("ring") and os.path.exists(m["ring"]) and os.path.exists(m["struct"])]
    if ms and len(sel) < n:
        r = dict(r); r["models"] = ms; sel.append(r)
print("rows in test:", [r["name"] for r in sel])
ids = [r["id"] for r in sel]
# subset of the frame: drop columns that are gaps in all selected rows (keeps the frame valid)
keep = [c for c in range(len(first[ids[0]])) if any(first[i][c] != "-" for i in ids)]
gap = {i: "".join(first[i][c] for c in keep) for i in ids}
# modules in the new columns
col2mod = {}
for name, a, b in mods:
    for c in range(a, b + 1): col2mod[c] = name
newmods = []
for k, c in enumerate(keep, 1):
    nm = col2mod[c + 1]
    if newmods and newmods[-1][0] == nm: newmods[-1][2] = k
    else: newmods.append([nm, k, k])
for r in sel:
    r["class"] = "test"; r["label"] = r["name"]; r["gene"] = r["name"].lower(); r["acc"] = r["accessions"][0][1]; r["header"] = r["name"]; r["info"] = []; r["flags"] = []
spec = {"pack": "cry_test", "system": "cry", "title": "test", "build": "test", "seg": "Cry toxin", "frames": [{"id": "struct", "label": "test frame", "kind": "sequence", "short": "x", "gapped": gap}], "order": ids, "rows": sel,
        "modules": [tuple(m) for m in newmods], "ref_row": ids[0], "default": {"auto": True, "row": ids[0], "why": "test"}, "node_text": "test", "inputs": {}}
t0 = time.time(); out = CP.build(spec); print("built in %.1fs" % (time.time() - t0)); print("strict columns", len(out["strict"]), "networks", out["nets"])
