"""Modules (domains) of the Cry Atlas.  Boundaries = consensus of the BPPRC domain annotation over the rows that have all three domains (D1, D2, D3), carried on the frame columns:
every residue of an annotated row votes for N-term / DI / DII / DIII / C-term (linker residues between two domains vote for the nearer domain); a monotone dynamic programme over the frame columns
(N-term <= DI <= DII <= DIII <= C-term) maximises the votes.  Every row then gets the module of the column of each of its residues (uniform boundaries, no per-row jitter).
Output work/modules.json: {frame, modules [[name, first col, last col]], per_row boundaries check}.     Run: py -I 15_modules.py [seq|struct]"""
import os, sys, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
frame = sys.argv[1] if len(sys.argv) > 1 else "seq"
rows = jload(os.path.join(WORK, "rows_seed.json")); R = {r["id"]: r for r in rows}
aln = dict(CC.read_fasta(os.path.join(WORK, "aln", "%s.fa" % frame)))
W = len(next(iter(aln.values())))
NAMES = ["N-term", "DI", "DII", "DIII", "C-term"]


def rng(s):
    if not s or s == "NOTFOUND": return None
    a, b = s.split("-"); return int(a), int(b)


votes = [[0] * 5 for _ in range(W + 1)]; nann = 0
for r in rows:
    d = [rng(r["bpprc_domains"].get(k)) for k in ("d1", "d2", "d3")]
    if any(x is None for x in d): continue
    nann += 1; (a1, b1), (a2, b2), (a3, b3) = d; cols = CC.col_index(aln[r["id"]])
    for p in range(1, len(r["seq"]) + 1):
        if p < a1: s = 0
        elif p <= b1: s = 1
        elif p < a2: s = 1 if p - b1 <= (a2 - b1) / 2.0 else 2
        elif p <= b2: s = 2
        elif p < a3: s = 2 if p - b2 <= (a3 - b2) / 2.0 else 3
        elif p <= b3: s = 3
        else: s = 4
        votes[cols[p - 1]][s] += 1
print("rows with all three BPPRC domains:", nann, "frame", frame, "width", W)
NEG = -10 ** 9
best = [[NEG] * 5 for _ in range(W + 1)]; back = [[0] * 5 for _ in range(W + 1)]
for s in range(5): best[0][s] = 0
for c in range(1, W + 1):
    for s in range(5):
        m, arg = max((best[c - 1][t], t) for t in range(s + 1))
        best[c][s] = m + votes[c][s]; back[c][s] = arg
s = max(range(5), key=lambda t: best[W][t]); path = [0] * (W + 1)
for c in range(W, 0, -1):
    path[c] = s; s = back[c][s]
mods = []
for c in range(1, W + 1):
    if mods and mods[-1][0] == path[c]: mods[-1][2] = c
    else: mods.append([path[c], c, c])
mods = [[NAMES[k], a, b] for k, a, b in mods]
print("modules:", mods)
# per-row check: where the module boundaries fall in the residues of annotated rows (mean offset to the BPPRC boundary)
off = collections.defaultdict(list)
col2mod = {}
for name, a, b in mods:
    for c in range(a, b + 1): col2mod[c] = name
for r in rows:
    d = [rng(r["bpprc_domains"].get(k)) for k in ("d1", "d2", "d3")]
    if any(x is None for x in d): continue
    cols = CC.col_index(aln[r["id"]]); mod = [col2mod[c] for c in cols]
    first = {n: next((i + 1 for i, m in enumerate(mod) if m == n), None) for n in NAMES}
    for k, n in enumerate(("DI", "DII", "DIII")):
        if first[n]: off[n].append(first[n] - d[k][0])
for n, v in off.items():
    v.sort(); print("  %-4s first residue of the module minus the BPPRC domain start: median %d, 10-90 %% %d..%d, n=%d" % (n, v[len(v) // 2], v[len(v) // 10], v[9 * len(v) // 10], len(v)))
cov = {n: sum(1 for r in rows if any(col2mod[c] == n for c in CC.col_index(aln[r["id"]]))) for n in NAMES}
print("rows with residues in each module:", cov)
jdump({"frame": frame, "width": W, "modules": mods, "n_annotated_rows": nann}, os.path.join(WORK, "modules_%s.json" % frame), indent=1)
