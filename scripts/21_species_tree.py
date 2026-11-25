"""GTDB species tree of the Cry Atlas (window 'Species tree' of the page) -> site/data/species_tree.json
Same file format as the species tree of the NODE-MSA Atlas (build_species_tree_20261105.py): the GTDB r232 bacterial tree pruned to the species of the Cry rows, one heat ring per Cry class,
one record per row (= distinct BPPRC sequence) in the cell of its species and class ring, linked to the row of the pack.
Species of a row = the organism of its GenPept record (first two words) matched to a GTDB species with a genome in the tree (letter suffixes of GTDB genera ignored: Bacillus_A thuringiensis = B. thuringiensis).
Not placed (listed in the window): no organism, 'Unknown', 'synthetic construct', genus-only names ('Bacillus sp.'), and Escherichia coli (the expression host named in the record, not the source of the gene), and organisms absent from GTDB.
GTDB files: GTDB_DIR (env) = folder with VERSION.txt, bac120.tree, bac120_taxonomy.tsv (about 140 MB; ~1 min).    Run: py -I 21_species_tree.py"""
import os, sys, re, json, collections, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
GT = os.environ.get("GTDB_DIR") or os.path.join(CRY, "gtdb_v232")
rows = jload(os.path.join(SITE, "data", "cry", "rows.json")); meta = jload(os.path.join(SITE, "data", "cry", "meta.json"))
classes = sorted({r["class"] for r in rows}, key=lambda c: -sum(1 for r in rows if r["class"] == c))
RINGS = classes
UNPLACED = {"", "unknown. unclassified.", "synthetic construct", "escherichia coli"}


def base_name(o):
    o = re.sub(r"[\[\]]", "", re.sub(r"^Candidatus\s+", "", (o or "").strip())); w = o.split()
    return " ".join(w[:2]) if len(w) >= 2 else o


def norm(n): return re.sub(r"_[A-Z]\b", "", n)


tree_txt = open(os.path.join(GT, "bac120.tree"), encoding="utf-8").read()
leaves = set(re.findall(r"[(,](RS_GCF_\d+\.\d|GB_GCA_\d+\.\d)(?=[:,)])", tree_txt))
leaf_sp, sp_count = {}, collections.Counter()
for l in open(os.path.join(GT, "bac120_taxonomy.tsv"), encoding="utf-8"):
    a, t = l.rstrip("\n").split("\t"); sp = t.split(";")[-1][3:]; sp_count[sp] += 1
    if a in leaves: leaf_sp.setdefault(sp, ("Bacteria", a, t))
norm_index = collections.defaultdict(set)
for sp in leaf_sp: norm_index[norm(sp)].add(sp)
print("GTDB leaves", len(leaves), "species with a leaf", len(leaf_sp), flush=True)


def gtdb_species(org):
    b = base_name(org)
    if b.lower() in UNPLACED or b.lower().endswith(" sp.") or " " not in b: return None
    cand = norm_index.get(norm(b)) or set()
    return min(cand, key=lambda c: (-sp_count[c], c)) if cand else None


SP = {}; unplaced = []; n_rec = 0; place_report = collections.Counter(); unplaced_orgs = collections.Counter()
for r in rows:
    sp = gtdb_species(r["species"]); names = (r.get("extra") or {}).get("names") or [r["label"]]
    if sp is None:
        unplaced.append({"pack": "cry", "id": r["id"], "org": r["species"] or "no organism in the record"}); unplaced_orgs[r["species"] or "(none)"] += 1; continue
    dom, leaf, lin = leaf_sp[sp]
    e = SP.setdefault(sp, dict(dom=dom, lin=lin, n=[0] * len(RINGS), ng=[0] * len(RINGS), cells={}, scanned=False)); k = RINGS.index(r["class"])
    prod = ", ".join(names[:4]) + (" +%d more with the same sequence" % (len(names) - 4) if len(names) > 4 else "")
    rec = [r["acc"], prod, "BPPRC entry (organism of its GenPept record)", "", 0, r["id"], len(names), "", [], len(names), names[:12]]
    e["cells"].setdefault(str(k), []).append(rec); e["n"][k] += 1; e["ng"][k] += len(names); n_rec += 1; place_report[sp] += 1
for e in SP.values():
    for k, recs in e["cells"].items(): recs.sort(key=lambda x: x[5])
print("placed rows", n_rec, "in", len(SP), "species:", dict(place_report.most_common()), "| not placed", len(unplaced), dict(unplaced_orgs.most_common(8)), flush=True)
labels = {sp: "%s (GTDB %s)" % (base_name(sp), sp) for sp in SP if norm(sp) != sp}


def parse_newick(s):
    parent, kids, name, blen = [-1], [[]], [""], [0.0]; cur = 0; i = 0; n = len(s)

    def new(p): parent.append(p); kids.append([]); name.append(""); blen.append(0.0); kids[p].append(len(parent) - 1); return len(parent) - 1
    while i < n:
        ch = s[i]
        if ch == "(": cur = new(cur); i += 1
        elif ch == ",": cur = new(parent[cur]); i += 1
        elif ch == ")": cur = parent[cur]; i += 1
        elif ch == ";": break
        elif ch == ":":
            j = i + 1
            while j < n and s[j] not in ",();": j += 1
            blen[cur] = float(s[i + 1:j]); i = j
        elif ch in " \n\r\t": i += 1
        else:
            if ch == "'": j = s.index("'", i + 1); name[cur] = s[i + 1:j]; i = j + 1
            else:
                j = i
                while j < n and s[j] not in ",():;": j += 1
                name[cur] = s[i:j]; i = j
    return parent, kids, name, blen


want = {leaf_sp[sp][1]: sp for sp in SP}
parent, kids, name, blen = parse_newick(tree_txt.strip())
keep = [False] * len(parent)
for i in range(len(parent)):
    if not kids[i] and name[i] in want: keep[i] = True
for i in range(len(parent) - 1, 0, -1):
    if keep[i]: keep[parent[i]] = True
out = {}; stack = [(0, False)]
while stack:
    v, done = stack.pop()
    if not done:
        stack.append((v, True))
        for c in kids[v]:
            if keep[c]: stack.append((c, False))
    else:
        ch = [out.pop(c) for c in kids[v] if keep[c]]
        if not kids[v]: out[v] = [want.get(name[v], name[v]), blen[v], []]
        elif len(ch) == 1: ch[0][1] += blen[v]; out[v] = ch[0]
        else: out[v] = ["", blen[v], ch]


def conv(x): return [round(x[1], 5), x[0]] if not x[2] else [round(x[1], 5), [conv(c) for c in x[2]]]


sys.setrecursionlimit(100000)
tree = conv(out[0]); print("tips", len(want), flush=True)
ver = open(os.path.join(GT, "VERSION.txt"), encoding="utf-8").read().replace("\n", " ").strip()
d = {"built": time.strftime("%Y-%m-%d"), "schema": 2, "record": ["protein id", "product", "basis", "assembly", "pack index", "row id", "n genomes", "IRM row id", "other Atlas rows [pack index, row id]", "n sequence variants", "protein ids"],
     "gtdb": ver, "rings": RINGS, "labels": labels, "packs": [{"id": "cry", "ring": "Cry"}], "trees": {"Bacteria": tree}, "species": SP, "unplaced": unplaced,
     "note": "One record per distinct BPPRC Cry sequence, in the species of the organism named in its GenPept record (not a genome scan): the species are those that BPPRC / NCBI name, so the tree is small and almost all entries are Bacillus thuringiensis. "
             "A ring is a Cry class of the pack; 'n genomes' counts the BPPRC names that have exactly this sequence."}
p = os.path.join(SITE, "data", "species_tree.json"); json.dump(d, open(p, "w", encoding="utf-8"), separators=(",", ":"), allow_nan=False)
print("written", p, os.path.getsize(p) // 1024, "KB")
