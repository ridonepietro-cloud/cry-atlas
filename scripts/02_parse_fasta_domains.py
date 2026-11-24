"""Read BPPRC FASTA downloads, locate each domain in the full-length sequence -> data/cry_seqs.tsv, data/cry_domains.tsv (1-based inclusive)"""
import os, re, csv, collections
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = os.path.join(ROOT, "raw", "fasta")
def rd(fn):
    d = collections.OrderedDict(); k = None
    for l in open(os.path.join(F, fn), encoding="utf8"):
        l = l.rstrip("\n")
        if l.startswith(">"): k = l[1:].strip(); d[k] = []
        elif k is not None: d[k].append(l.strip())
    return collections.OrderedDict((k, "".join(v)) for k, v in d.items())
full = rd("cry_full.fasta"); hol = rd("holotype.fasta"); mut = rd("mutant.fasta")
dom = {i: rd("cry_domain%d.fasta" % i) for i in (1, 2, 3)}
print("full", len(full), "holotype", len(hol), "mutant", list(mut))
tab = {r["name"]: r for r in csv.DictReader(open(os.path.join(ROOT, "data", "bpprc_cry_table.tsv"), encoding="utf8"), delimiter="\t")}
missing = [n for n in tab if n not in full]; extra = [n for n in full if n not in tab]
print("table rows without sequence:", missing); print("sequence without table row:", extra)
# chars
bad = {n: set(s) - set("ACDEFGHIKLMNPQRSTVWY") for n, s in full.items()}
print("non-standard residues:", {n: sorted(v) for n, v in bad.items() if v})
import statistics
L = [len(s) for s in full.values()]
print("length min/median/max", min(L), statistics.median(L), max(L))
print("lengths <300:", [(n, len(s)) for n, s in full.items() if len(s) < 300][:40])
print("lengths >1300:", [(n, len(s)) for n, s in full.items() if len(s) > 1300][:40])
# domain location
rows = []; nfail = collections.Counter(); multi = 0
for n, s in full.items():
    r = {"name": n, "len": len(s)}
    for i in (1, 2, 3):
        ds = dom[i].get(n + "_d%d" % i)
        if ds is None: r["d%d" % i] = ""; continue
        p = s.find(ds)
        if p < 0: nfail[i] += 1; r["d%d" % i] = "NOTFOUND"; continue
        if s.count(ds) > 1: multi += 1
        r["d%d" % i] = "%d-%d" % (p + 1, p + len(ds))
    rows.append(r)
print("domain-not-found:", dict(nfail), "multi-occurrence:", multi)
nod = collections.Counter(tuple(bool(r["d%d" % i]) for i in (1, 2, 3)) for r in rows); print("domain presence patterns:", nod)
with open(os.path.join(ROOT, "data", "cry_domains_bpprc.tsv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["name", "len", "d1", "d2", "d3"], delimiter="\t"); w.writeheader(); w.writerows(rows)
# boundary statistics for Cry1Aa1
for r in rows[:1]: print(r)
ex = {r["name"]: r for r in rows}
for n in ("Cry1Aa1", "Cry3Aa1", "Cry2Aa1", "Cry4Aa1", "Cry11Aa1", "Cry8Ea1", "Cry5Ba1"): print(ex.get(n))
print("holotype in full?", sum(1 for k in hol if k in full), "of", len(hol)); print(list(hol)[:12])
