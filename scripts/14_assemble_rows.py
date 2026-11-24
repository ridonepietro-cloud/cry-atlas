"""Assemble the rows of the Cry Atlas pack (work/rows.json): metadata, structure models, specificity.
Models per row (any of them may be missing):
  experimental : PDB entities >= 97 % identical to the row over >= 70 % of the entity (work/pdb/matches.json), first chain with residues; scheme = SEQRES scheme of the mmCIF
  AF2 (AlphaFold DB, exact sequence)            : work/afdb/map.json, RAW/afdb
  AF2 donor (flagged)                           : rows without any of the above get the AlphaFold DB model of the most similar row with a model when it is >= 97 % identical over >= 95 % of the row
Primary model of a row = the one that places most residues of the row sequence (experimental wins a tie): full-length models first (the Atlas rule).
Output work/rows.json, work/models_report.tsv.    Run: py -I 14_assemble_rows.py"""
import os, sys, re, json, csv, collections, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
rows = jload(os.path.join(WORK, "rows_seed.json")); R = {r["id"]: r for r in rows}
mp = jload(os.path.join(WORK, "afdb", "map.json")); matches = jload(os.path.join(WORK, "pdb", "matches.json"), {}); spec = jload(os.path.join(WORK, "specificity.json"), {})
import hashlib
afdb_of = {}
for md5, v in mp.items():
    if v.get("afdb"): afdb_of[md5] = v
RING = os.path.join(RAW, "ring_out")


def afdb_file(a): return os.path.join(RAW, "afdb", os.path.basename(a["pdbUrl"]))


from Bio.Align import PairwiseAligner
al = PairwiseAligner(); al.mode = "global"; al.match_score = 2; al.mismatch_score = -1; al.open_gap_score = -6; al.extend_gap_score = -1; al.end_gap_score = 0


def ident_cov(a, b):
    x = al.align(a, b)[0]; same = cov = 0
    for (r0, r1), (c0, c1) in zip(*x.aligned):
        cov += r1 - r0; same += sum(1 for q in range(r1 - r0) if a[r0 + q] == b[c0 + q])
    return (same / cov if cov else 0.0), cov


def kset(s, k=5): return {s[i:i + k] for i in range(len(s) - k + 1)}


for r in rows:
    r["models"] = []; md5 = hashlib.md5(r["seq"].encode()).hexdigest().upper()
    a = afdb_of.get(md5)
    r["uniprot"] = mp.get(md5, {}).get("uniprot", [])
    if a:
        f = afdb_file(a["afdb"])
        r["models"].append({"label": "AlphaFold DB %s (v%s, mean pLDDT %s)" % (a["afdb"]["entry"], a["afdb"].get("version"), a["afdb"].get("plddt")), "kind": "AF2", "struct": f, "chains": ["A"],
                            "ring": os.path.join(RING, os.path.basename(f) + "_ringEdges"), "pdb": None, "uniprot": a["afdb"]["acc"]})
# experimental
for pid, ents in sorted(matches.items()):
    cif = os.path.join(RAW, "pdb", pid + ".cif")
    if not os.path.exists(cif): continue
    sch = None
    for en in ents:
        for rid, ident, frac in en["rows"]:
            if ident < 0.97: continue
            r = R[rid]
            if sch is None: sch = CC.cif_scheme(cif)
            chains = [c for c in en["chains"]]
            r["models"].append({"label": "Experimental %s (entity %s, %d aa%s)" % (pid, en["entity"], en["len"], "" if ident >= 0.9999 else ", %.1f %% identical to this entry" % (100 * ident)), "kind": "experimental", "struct": cif, "chains": chains, "pdb": pid,
                                "ring": os.path.join(RING, pid + ".cif_ringEdges"), "scheme_ok": True, "ident": ident, "entity_len": en["len"], "title": en["title"]})
# donors
have = [r for r in rows if r["models"]]; ks = {r["id"]: kset(r["seq"]) for r in have}
n_donor = 0
for r in rows:
    if r["models"]: continue
    k = kset(r["seq"]); cand = sorted(((len(k & ks[d["id"]]) / max(1, min(len(k), len(ks[d["id"]]))), d["id"]) for d in have if any(m["kind"] == "AF2" for m in d["models"])), reverse=True)[:4]
    best = None
    for sh, did in cand:
        if sh < 0.5: break
        ident, cov = ident_cov(r["seq"], R[did]["seq"])
        if ident >= 0.97 and cov >= 0.95 * len(r["seq"]) and (best is None or ident > best[0]): best = (ident, cov, did)
    if best:
        d = R[best[2]]; m = next(m for m in d["models"] if m["kind"] == "AF2"); mm = dict(m); mm["donor"] = d["id"]; mm["donor_name"] = d["name"]; mm["donor_ident"] = best[0]
        mm["label"] = "AlphaFold DB model of %s (%.1f %% identical to this entry; not a model of this sequence)" % (d["name"], 100 * best[0])
        r["models"].append(mm); n_donor += 1
print("rows", len(rows), "with AFDB exact", sum(1 for r in rows if any(m["kind"] == "AF2" and "donor" not in m for m in r["models"])),
      "with experimental", sum(1 for r in rows if any(m["kind"] == "experimental" for m in r["models"])), "donor-only", n_donor, "no model", sum(1 for r in rows if not r["models"]))
# specificity layer
for r in rows:
    t = {"yes": collections.Counter(), "no": collections.Counter(), "sp_yes": set(), "sp_no": set()}; got = False
    for n in r["names"]:
        s = spec.get(n)
        if not s: continue
        got = True
        for o, c in s["yes"].items(): t["yes"][o] += c
        for o, c in s["no"].items(): t["no"][o] += c
        t["sp_yes"] |= set(s["sp_yes"]); t["sp_no"] |= set(s["sp_no"])
    r["target"] = {"yes": dict(t["yes"]), "no": dict(t["no"]), "sp_yes": sorted(t["sp_yes"]), "sp_no": sorted(t["sp_no"])} if got else None
jdump(rows, os.path.join(WORK, "rows.json"), indent=0)
with open(os.path.join(WORK, "models_report.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["id", "name", "len", "n_models", "models"])
    for r in rows: w.writerow([r["id"], r["name"], r["len"], len(r["models"]), " | ".join(m["label"] for m in r["models"])])
