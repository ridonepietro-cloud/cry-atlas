"""Activity breadth and clades (2026-11-27): (a) how broad is the killing of nematode-active vs insect-active Cry toxins in the BPPRC specificity data; (b) do the nematode-active toxins form a clade in the Atlas trees.
Input : raw/specificity/association.csv (BPPRC export), site/data/cry/rows.json + trees.json.   Output: printed tables + work/specificity_breadth.json.     Run: py -I 22_specificity_breadth_trees.py"""
import os, sys, csv, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cry_common import *
rows = jload(os.path.join(SITE, "data", "cry", "rows.json")); name2row = {n: r["id"] for r in rows for n in (r.get("extra") or {}).get("names", [])}
rd = [x for x in csv.DictReader(open(os.path.join(CRY, "raw", "specificity", "association.csv"), encoding="utf-8", newline="")) if x["name"] in name2row]
ORDERS = ["Rhabditida", "Lepidoptera", "Coleoptera", "Diptera"]
# per toxin (BPPRC name) and order: species with an active record, species with an inactive record
act = collections.defaultdict(lambda: collections.defaultdict(lambda: {"yes": set(), "no": set(), "rec_yes": 0, "rec_no": 0}))
for x in rd:
    a = (x["activity"] or "").strip().lower(); o = x["target_order"].strip(); sp = x["target_species"].strip()
    if a not in ("yes", "no") or not o: continue
    d = act[x["name"]][o]; d[a].add(sp); d["rec_" + a] += 1
res = {}
print("order        toxins tested | toxins with >=1 active species | of those: also >=1 inactive species | species-toxin pairs: active / inactive (both = counted in each)")
for o in ORDERS:
    tox = [n for n in act if o in act[n]]; pos = [n for n in tox if act[n][o]["yes"]]; mixed = [n for n in pos if act[n][o]["no"]]
    only_neg = [n for n in tox if not act[n][o]["yes"] and act[n][o]["no"]]
    ya = sum(len(act[n][o]["yes"]) for n in tox); na = sum(len(act[n][o]["no"]) for n in tox)
    res[o] = {"toxins_tested": len(tox), "toxins_active": len(pos), "active_with_inactive_species": len(mixed), "only_negative": len(only_neg), "pairs_active": ya, "pairs_inactive": na,
              "mixed_names": sorted(mixed)[:40], "active_names": sorted(pos)}
    print("%-12s %5d        | %5d                          | %5d (%.0f %%)                        | %d / %d   (toxins with only negative records: %d)" % (o, len(tox), len(pos), len(mixed), 100.0 * len(mixed) / max(1, len(pos)), ya, na, len(only_neg)))
print()
nem = res["Rhabditida"]["active_names"]
print("nematode-active BPPRC names (%d):" % len(nem))
for n in nem:
    d = act[n]["Rhabditida"]; ins = {o: len(act[n][o]["yes"]) for o in act[n] if o != "Rhabditida" and act[n][o]["yes"]}
    print("  %-10s nematode species active %s | inactive %s | also active in insect orders %s" % (n, sorted(d["yes"]), sorted(d["no"]), ins or "none recorded"))
# toxins active on nematodes AND on insects
both = [n for n in nem if any(o != "Rhabditida" and act[n][o]["yes"] for o in act[n])]
print("nematode-active names that also have an active insect record:", both)
nem_rows = sorted({name2row[n] for n in nem}); json.dump({"orders": res, "nematode_rows": nem_rows, "nematode_names": nem}, open(os.path.join(WORK, "specificity_breadth.json"), "w"), indent=0)
# insects: how patchy? distribution of fraction of tested species that are active, per toxin (toxins with >= 3 tested species in the order)
for o in ORDERS:
    fr = []
    for n in act:
        d = act[n].get(o)
        if not d: continue
        tested = d["yes"] | d["no"]
        if len(tested) >= 3: fr.append(len(d["yes"]) / float(len(tested)))
    if fr:
        fr.sort(); print("%-12s toxins with >= 3 tested species: %d | median fraction of tested species killed %.2f | toxins that killed ALL tested: %d | killed NONE or some but not all: %d" % (o, len(fr), fr[len(fr) // 2], sum(1 for f in fr if f == 1.0), sum(1 for f in fr if f < 1.0)))

# ---- (b) clade test on the trees (tree code of the ASR tools: Newick parser + midpoint root; ete3 does not run on this Python)
sys.path.insert(0, os.path.join(HERE, "asr"))
import asr_lib as L
T = jload(os.path.join(SITE, "data", "cry", "trees.json"))["trees"]
tgt = set(nem_rows); print(chr(10) + "nematode-active rows:", len(tgt), sorted(tgt))
cls = {r["id"]: r["class"] for r in rows}
for key in ("struct_A", "seq_A", "struct_DI", "struct_DII", "struct_DIII", "seq_DI", "seq_DII", "seq_DIII"):
    root = L.midpoint_root(L.parse_newick(T[key]["newick"])); leaves_all = {l.name for l in L.leaves(root)}; tt = tgt & leaves_all
    if len(tt) < 2: continue
    lv = {}
    for nd in L.postorder(root): lv[id(nd)] = {nd.name} if not nd.kids else set().union(*[lv[id(c)] for c in nd.kids])
    cand = [nd for nd in L.postorder(root) if tt <= lv[id(nd)]]; lca = min(cand, key=lambda nd: len(lv[id(nd)])); inside = lv[id(lca)]
    best = max((len(lv[id(nd)]) for nd in L.postorder(root) if lv[id(nd)] <= tt), default=0)
    classes = collections.Counter(cls[i] for i in inside)
    print("%-12s leaves %d | nematode-active leaves %d | smallest clade holding all of them: %d leaves (%.0f %% nematode-active), SH-aLRT %s | largest clade made only of nematode-active leaves: %d | classes inside: %s" % (key, len(leaves_all), len(tt), len(inside), 100.0 * len(tt) / len(inside), lca.label or "-", best, dict(classes.most_common(4))))

# ---- (c) the nematode families only: Cry5 / Cry12 / Cry13 / Cry14 / Cry21 / Cry58 / Cry95 (class 'Cry21 / Cry5 / Cry14 +3' plus Cry95) against the rest
fam_of = {r["id"]: (r.get("extra") or {}).get("family") for r in rows}; NEMF = {"Cry5", "Cry12", "Cry13", "Cry14", "Cry21", "Cry58", "Cry95"}
core_all = {i for i, f in fam_of.items() if f in NEMF}; core_act = tgt & core_all; odd = tgt - core_all
print(chr(10) + "rows of the families %s: %d | nematode-active among them %d | nematode-active rows outside these families: %s" % (sorted(NEMF), len(core_all), len(core_act), sorted((i, fam_of[i]) for i in odd)))
print("rows of these families WITHOUT a nematode-active record: %s" % sorted((i, fam_of[i]) for i in core_all - core_act))
for key in ("struct_A", "seq_A", "struct_DI", "struct_DII", "struct_DIII"):
    root = L.midpoint_root(L.parse_newick(T[key]["newick"])); lv = {}
    for nd in L.postorder(root): lv[id(nd)] = {nd.name} if not nd.kids else set().union(*[lv[id(c)] for c in nd.kids])
    leaves_all = lv[id(root)]; fam = core_all & leaves_all
    cand = [nd for nd in L.postorder(root) if fam <= lv[id(nd)]]; lca = min(cand, key=lambda nd: len(lv[id(nd)])); inside = lv[id(lca)]
    print("%-12s family rows in tree %d: smallest clade holding all of them has %d leaves (%.0f %% are rows of these families), SH-aLRT %s | outsiders inside: %s" % (key, len(fam), len(inside), 100.0 * len(fam) / len(inside), lca.label or "-", sorted((i, fam_of[i]) for i in inside - fam)[:8]))
