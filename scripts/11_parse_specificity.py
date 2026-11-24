"""BPPRC specificity export (raw/specificity/association.csv) -> work/specificity.json: per Cry name, the target orders / species with activity Yes / No, LC50 records, sources.
The curator column (data_entered_by) is not carried over.   Run: py -I 11_parse_specificity.py"""
import os, sys, csv, json, re, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
rows = jload(os.path.join(WORK, "rows_seed.json"))
names = {n for r in rows for n in r["names"]}
rd = list(csv.DictReader(open(os.path.join(CRY, "raw", "specificity", "association.csv"), encoding="utf-8", newline="")))
print("association rows", len(rd), "columns", list(rd[0].keys()))
cry = [x for x in rd if x["name"] in names]
print("rows for Cry names in the table:", len(cry), "distinct names", len({x["name"] for x in cry}), "of", len(names))
print("activity values:", collections.Counter(x["activity"] for x in cry))
print("target orders:", collections.Counter(x["target_order"] for x in cry).most_common(20))
out = {}
for x in cry:
    d = out.setdefault(x["name"], {"yes": collections.Counter(), "no": collections.Counter(), "sp_yes": set(), "sp_no": set(), "partner": set(), "lc50": [], "src": set()})
    act = (x["activity"] or "").strip().lower()
    o = x["target_order"].strip() or "unknown"
    k = "yes" if act == "yes" else "no" if act == "no" else None
    if k is None: continue
    d[k][o] += 1; d["sp_" + k].add(x["target_species"].strip())
    if x["partnerprotein"].strip().lower() == "yes": d["partner"].add(x["partnerprotein_textbox"].strip())
    if x["lc50"].strip(): d["lc50"].append([x["target_species"].strip(), x["lc50"].strip(), x["units"].strip()])
    src = x["publication"].strip() or x["other_citations"].strip()
    if src: d["src"].add(src)
ser = {}
for n, d in out.items():
    ser[n] = {"yes": dict(d["yes"]), "no": dict(d["no"]), "sp_yes": sorted(d["sp_yes"]), "sp_no": sorted(d["sp_no"]), "partner": sorted(p for p in d["partner"] if p and p != "Not applicable"), "lc50": d["lc50"][:12], "src": sorted(d["src"])[:6]}
jdump(ser, os.path.join(WORK, "specificity.json"), indent=0)
act = collections.Counter(o for d in ser.values() for o in d["yes"])
print("names with an active record:", sum(1 for d in ser.values() if d["yes"]), "| active-order counts (names):", collections.Counter(o for d in ser.values() for o in d["yes"]).most_common(12))
print("rows (distinct sequences) with at least one name in the table:", sum(1 for r in rows if any(n in ser for n in r["names"])), "of", len(rows))
