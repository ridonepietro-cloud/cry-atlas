"""Rows of the Cry Atlas: one row per DISTINCT BPPRC sequence (exact duplicates are merged and all their names listed).
Inputs: data/bpprc_cry_table.tsv, data/cry_seq_clean.json, data/cry_domains_bpprc.tsv, raw/fasta/holotype.fasta, work/ncbi/genpept.json, work/ncbi/compare.tsv
Output: work/rows_seed.json  (list of row dicts, ids C000.. in natural BPPRC name order)    Run: py -I 09_build_rows.py"""
import os, sys, re, csv, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *

tab = list(csv.DictReader(open(os.path.join(DATA, "bpprc_cry_table.tsv"), encoding="utf8"), delimiter="\t"))
T = {r["name"]: r for r in tab}
seqs = jload(os.path.join(DATA, "cry_seq_clean.json")); issues = seqs["issues"]; seqs = seqs["seqs"]
gp = jload(os.path.join(WORK, "ncbi", "genpept.json"), {})
dom = {r["name"]: r for r in csv.DictReader(open(os.path.join(DATA, "cry_domains_bpprc.tsv"), encoding="utf8"), delimiter="\t")}
hol = set(l[1:].strip() for l in open(os.path.join(CRY, "raw", "fasta", "holotype.fasta"), encoding="utf8") if l.startswith(">"))
cmp_ = {r["name"]: r for r in csv.DictReader(open(os.path.join(WORK, "ncbi", "compare.tsv"), encoding="utf8"), delimiter="\t")}
# the second record of Cry3Aa10 is the one equal to the BPPRC detail page / NCBI length (652): choose by NCBI length when a name has two different records
raw_recs = collections.defaultdict(list); k = None
for l in open(os.path.join(CRY, "raw", "fasta", "cry_full.fasta"), encoding="utf8"):
    l = l.rstrip("\n")
    if l.startswith(">"): k = l[1:].strip(); raw_recs[k].append([])
    elif k is not None: raw_recs[k][-1].append(l)
for n, lst in raw_recs.items():
    ss = {re.sub(r"\s+", "", "".join(x)).upper() for x in lst}
    if len(ss) > 1:
        g = gp.get(T[n]["accession"]) or {}
        want = g.get("length")
        pick = [s for s in ss if s == g.get("seq")] or [s for s in ss if len(s) == want]
        print("name with two different records:", n, sorted(len(s) for s in ss), "NCBI length", want, "-> picked", [len(s) for s in pick])
        if pick: seqs[n] = pick[0]; issues.setdefault(n, []).append("two different records in the BPPRC download (%s aa); the one with the NCBI length is used" % "/".join(str(len(s)) for s in sorted(ss)))


def nkey(name):
    m = re.match(r"Cry(\d+)([A-Za-z]*)(\d*)(.*)", name)
    if not m: return (999, name, 0, "")
    return (int(m.group(1)), m.group(2).lower() if m.group(2) else "", int(m.group(3) or 0), m.group(4))


def family(name):
    m = re.match(r"(Cry\d+)", name); return m.group(1) if m else "other"


by = collections.OrderedDict()
for n in sorted(seqs, key=nkey):
    by.setdefault(seqs[n], []).append(n)
rows = []
for i, (s, names) in enumerate(by.items()):
    n0 = names[0]; t0 = T[n0]
    orgs = collections.Counter(); strains = []; plasmids = []; defs = []; acc = []; coded = []
    for n in names:
        a = T[n]["accession"]; g = gp.get(a); acc.append([n, a])
        if g and not g.get("missing"):
            if g.get("organism"): orgs[g["organism"]] += 1
            for st in g.get("strain", []): strains.append(st)
            for p in g.get("plasmid", []): plasmids.append(p)
            defs.append(g.get("definition", ""))
            if g.get("coded_by"): coded.append([n, g["coded_by"]])
    years = [int(T[n]["year"]) for n in names if T[n]["year"].isdigit()]
    d = dom.get(n0, {})
    flags, info = [], []
    for n in names:
        for x in issues.get(n, []): info.append("%s: %s" % (n, x))
        c = cmp_.get(n)
        if c and c["result"] not in ("identical", "no NCBI protein record"):
            info.append("%s: NCBI protein %s differs from the BPPRC sequence (%s)" % (n, T[n]["accession"], c["detail"] or c["result"]))
    if len(names) > 1: info.append("%d BPPRC entries have exactly this sequence: %s" % (len(names), ", ".join(names)))
    r = {"id": "C%03d" % i, "name": n0, "names": names, "family": family(n0), "len": len(s), "seq": s, "accessions": acc,
         "species": (orgs.most_common(1)[0][0] if orgs else ""), "strains": sorted(set(strains)), "plasmids": sorted(set(plasmids)), "definitions": sorted(set(defs))[:3],
         "year": min(years) if years else None, "holotype": any(n in hol for n in names), "old_names": sorted({T[n]["old_name"] for n in names if T[n]["old_name"] and T[n]["old_name"] != n}),
         "other_names": sorted({T[n]["other_name"] for n in names if T[n]["other_name"]}),
         "bpprc_domains": {k: d.get(k, "") for k in ("d1", "d2", "d3")}, "coded_by": coded[:3], "flags": flags, "info": info}
    rows.append(r)
print("rows", len(rows), "from", len(seqs), "names")
print("families", len({r["family"] for r in rows}), collections.Counter(r["family"] for r in rows).most_common(12))
print("rows with organism", sum(1 for r in rows if r["species"]), "with holotype", sum(1 for r in rows if r["holotype"]))
print("top species", collections.Counter(r["species"] for r in rows).most_common(12))
jdump(rows, os.path.join(WORK, "rows_seed.json"), indent=0)
# no-sequence entries
nos = [r for r in tab if r["name"] not in seqs]
jdump(nos, os.path.join(DATA, "no_sequence_entries.json"), indent=0)
print("entries with no sequence in the BPPRC download:", [(r["name"], r["accession"]) for r in nos])
