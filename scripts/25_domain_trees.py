"""Domain trees compared entry by entry (domain swaps / recombination between Domain I, II and III).
For every entry x and every pair of domain trees: the K = 10 nearest relatives of x (patristic distance) in tree 1 and in tree 2; discordance d = 1 - Jaccard(neighbourhood 1, neighbourhood 2).
Calibration of the noise: the SAME domain from the structure-based and the sequence-based alignment gives two independent trees of the same history; their discordance d_rep (per entry) is the background.
Threshold = 90th percentile of d_rep over all entries and domains.  Category of an entry (structure-based trees): 'concordant' = no pair above the threshold; 'Domain X differs' = the two pairs that contain X are above, the other pair is not;
'all three differ' = all pairs above; 'not in all trees' = the entry is missing from a domain tree (too few residues of a domain).
Also the tree-to-tree normalised Robinson-Foulds distance (unrooted splits on the common entries) for DI-DII, DI-DIII, DII-DIII and for the replicates.
Output work/domain_discordance.json, docs/domain_tree_discordance.md, docs/domain_tree_discordance.tsv.      Run: py -I 25_domain_trees.py"""
import os, sys, json, csv, collections
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "asr"))
from cry_common import *
import asr_lib as L
K = 10
rows = jload(os.path.join(SITE, "data", "cry", "rows.json")); lab = {r["id"]: r["label"].split(" ")[0] for r in rows}; fam = {r["id"]: (r.get("extra") or {}).get("family") for r in rows}
cls = {r["id"]: r["class"] for r in rows}; T = jload(os.path.join(SITE, "data", "cry", "trees.json"))["trees"]


def prep(key):
    root = L.midpoint_root(L.parse_newick(T[key]["newick"])); nodes = L.postorder(root); leaves = [n for n in nodes if not n.kids]; names = [n.name for n in leaves]; pos = {n: i for i, n in enumerate(names)}
    nl = len(names); D = np.zeros((nl, nl), dtype=np.float32); lset = {}; up = {}
    for nd in nodes:
        if not nd.kids: lset[id(nd)] = [pos[nd.name]]; up[id(nd)] = np.array([nd.len], dtype=np.float32); continue
        kids = nd.kids; parts = [(lset[id(c)], up[id(c)]) for c in kids]
        for a in range(len(parts)):
            for b in range(a + 1, len(parts)):
                (la, ua), (lb, ub) = parts[a], parts[b]; D[np.ix_(la, lb)] = ua[:, None] + ub[None, :]; D[np.ix_(lb, la)] = (ua[:, None] + ub[None, :]).T
        allL = [i for l, u in parts for i in l]; allU = np.concatenate([u for l, u in parts]) + np.float32(nd.len); lset[id(nd)] = allL; up[id(nd)] = allU
    return names, D, root


def neighbours(names, D, common, k=K):
    idx = [i for i, n in enumerate(names) if n in common]; sub = D[np.ix_(idx, idx)].copy(); np.fill_diagonal(sub, np.inf); nm = [names[i] for i in idx]
    nn = {}
    for r in range(len(idx)):
        order = np.argsort(sub[r], kind="stable")[:k]; nn[nm[r]] = [nm[j] for j in order]
    return nn


def jacc(a, b): return len(set(a) & set(b)) / float(len(set(a) | set(b)))


def splits(root, common):
    full = 0; bit = {}
    for i, n in enumerate(sorted(common)): bit[n] = 1 << i; full |= 1 << i
    S = set(); mask = {}
    for nd in L.postorder(root):
        if not nd.kids: mask[id(nd)] = bit.get(nd.name, 0); continue
        m = 0
        for c in nd.kids: m |= mask[id(c)]
        mask[id(nd)] = m
        if m and m != full and bin(m).count("1") > 1 and bin(full ^ m).count("1") > 1:
            S.add(m if not (m & 1) else full ^ m)
    return S


keys = {"DI": ("struct_DI", "seq_DI"), "DII": ("struct_DII", "seq_DII"), "DIII": ("struct_DIII", "seq_DIII")}
P = {k: prep(k) for ks in keys.values() for k in ks}
allsets = [set(P[k][0]) for k in P]; common_all = set.intersection(*allsets); print("entries in all six domain trees:", len(common_all), "| in struct trees:", len(set.intersection(*[set(P[keys[d][0]][0]) for d in keys])))
NN = {k: neighbours(P[k][0], P[k][1], common_all) for k in P}
import re
G2 = lambda nm: (re.match(r"(Cry\d+[A-Z])", nm) or [None, nm])[1]                 # secondary-rank group of a BPPRC name (Cry1A, Cry8B ...): only a yardstick of similarity for the neighbours, not used to build any tree
def majority(tree, x):
    g = collections.Counter(G2(lab[n]) for n in NN[tree][x]); top, n = g.most_common(1)[0]; return top, n / float(K)
MINSHARE = 0.6
# noise: the same domain from the structure-based and the sequence-based alignment (two trees of the same history)
rep_diff = {}
for d in keys:
    ok = diff = 0
    for x in common_all:
        (g1, s1), (g2, s2) = majority(keys[d][0], x), majority(keys[d][1], x)
        if s1 >= MINSHARE and s2 >= MINSHARE: ok += 1; diff += g1 != g2
    rep_diff[d] = (diff, ok); print("replicate trees of %s: entries with a clear neighbourhood in both %d, majority group differs in %d (%.1f %%)" % (d, ok, diff, 100.0 * diff / max(1, ok)))
cat = {}; nn1 = {}; maj = {}; score = {}
pairs = [("DI", "DII"), ("DI", "DIII"), ("DII", "DIII")]
for x in common_all:
    m = {d: majority(keys[d][0], x) for d in keys}; maj[x] = m; nn1[x] = {d: NN[keys[d][0]][x][0] for d in keys}
    if min(v[1] for v in m.values()) < MINSHARE: cat[x] = "mixed neighbourhood"; score[x] = 0.0; continue
    g = {d: m[d][0] for d in keys}; vals = list(g.values()); nd = len(set(vals))
    if nd == 1: cat[x] = "concordant"
    elif nd == 3: cat[x] = "all three differ"
    else:
        odd = [d for d in keys if vals.count(g[d]) == 1][0]; cat[x] = {"DI": "Domain I differs", "DII": "Domain II differs", "DIII": "Domain III differs"}[odd]
    score[x] = (nd - 1) / 2.0
disc = {p: {x: 1 - jacc(NN[keys[p[0]][0]][x], NN[keys[p[1]][0]][x]) for x in common_all} for p in pairs}
pool = np.array([1 - jacc(NN[keys[d][0]][x], NN[keys[d][1]][x]) for d in keys for x in common_all]); thr = float(np.percentile(pool, 90))
cc = collections.Counter(cat.values()); print("categories (majority secondary-rank group of the 10 nearest relatives, clear neighbourhood >= %.0f %%):" % (100 * MINSHARE), dict(cc))
assess = sum(v for k, v in cc.items() if k != "mixed neighbourhood"); print("assessable entries %d: concordant %.0f %%, one domain differs %.0f %%, all three differ %.0f %%" % (assess, 100.0 * cc["concordant"] / assess, 100.0 * (cc["Domain I differs"] + cc["Domain II differs"] + cc["Domain III differs"]) / assess, 100.0 * cc["all three differ"] / assess))
for p in pairs:
    ok = diff = 0
    for x in common_all:
        a_, b_ = maj[x][p[0]], maj[x][p[1]]
        if a_[1] >= MINSHARE and b_[1] >= MINSHARE: ok += 1; diff += a_[0] != b_[0]
    print("domain pair %s-%s: majority group differs in %d of %d entries with a clear neighbourhood (%.1f %%; replicate noise %s %.1f %% / %s %.1f %%)" % (p[0], p[1], diff, ok, 100.0 * diff / max(1, ok), p[0], 100.0 * rep_diff[p[0]][0] / max(1, rep_diff[p[0]][1]), p[1], 100.0 * rep_diff[p[1]][0] / max(1, rep_diff[p[1]][1])))
# normalised RF between the structure-based domain trees and the replicates
spl = {k: splits(P[k][2], common_all) for k in P}
def nrf(a, b): return len(spl[a] ^ spl[b]) / float(len(spl[a]) + len(spl[b]))
rf = {"%s-%s" % p: nrf(keys[p[0]][0], keys[p[1]][0]) for p in pairs}; rf_rep = {d: nrf(*keys[d]) for d in keys}
print("normalised RF between domain trees (structure frame):", {k: round(v, 2) for k, v in rf.items()}, "| replicates (structure vs sequence frame, same domain):", {k: round(v, 2) for k, v in rf_rep.items()})
# entries with a domain that differs, strongest first (share of the odd group's neighbours, closest relative per domain)
top = [x for x in sorted(common_all, key=lambda x: (-score[x], -min(v[1] for v in maj[x].values()))) if cat[x] not in ("concordant", "mixed neighbourhood")][:30]
by_cls = collections.defaultdict(lambda: collections.Counter())
for x in common_all: by_cls[cls[x]][cat[x]] += 1
os.makedirs(os.path.join(CRY, "docs"), exist_ok=True)
with open(os.path.join(CRY, "docs", "domain_tree_discordance.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["row", "entry", "family", "class", "category", "group_DI", "group_DII", "group_DIII", "share_DI", "share_DII", "share_DIII", "closest_DI", "closest_DII", "closest_DIII"])
    for x in sorted(common_all, key=lambda x: (-score[x], lab[x])):
        w.writerow([x, lab[x], fam[x], cls[x], cat[x]] + [maj[x][d][0] for d in ("DI", "DII", "DIII")] + [round(maj[x][d][1], 2) for d in ("DI", "DII", "DIII")] + [lab[nn1[x][d]] for d in ("DI", "DII", "DIII")])
json.dump({"K": K, "min_share": MINSHARE, "replicate_noise": {d: [rep_diff[d][0], rep_diff[d][1]] for d in keys}, "category": cat, "closest": nn1, "groups": {x: {d: [maj[x][d][0], round(maj[x][d][1], 2)] for d in keys} for x in common_all},
           "rf": rf, "rf_replicate": rf_rep, "counts": dict(cc), "n_entries": len(common_all)}, open(os.path.join(WORK, "domain_discordance.json"), "w"), indent=0)
print("most discordant entries (mean discordance; closest relative in the DI / DII / DIII tree):")
for x in top[:25]: print("  %-12s %-20s | group of the nearest relatives DI / DII / DIII: %s / %s / %s | closest: %s / %s / %s" % (lab[x], cat[x], maj[x]["DI"][0], maj[x]["DII"][0], maj[x]["DIII"][0], lab[nn1[x]["DI"]], lab[nn1[x]["DII"]], lab[nn1[x]["DIII"]]))
print("category by class:");
for c, v in by_cls.items(): print("  %-28s %s" % (c, dict(v)))
