"""Class of every row (colour / filter group of the page, at most 8): Cry primary ranks (Cry1, Cry2, ...) grouped by the similarity of their toxin cores.
Distance between two rows = p-distance over the columns of DI + DII + DIII where both have residues (segment frame work/aln/seqseg.fa); distance between two families = mean over row pairs;
average-linkage clustering of the families, cut into N_CLASS groups; a group is named by its member families with the most rows ("Cry1 / Cry9 ...").  Output work/classes.json {family: class}.
Run: py -I 17_classes.py [frame]"""
import os, sys, json, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
frame = sys.argv[1] if len(sys.argv) > 1 else "seqseg"; N_CLASS = int(os.environ.get("CRY_NCLASS", "8"))
rows = jload(os.path.join(WORK, "rows_seed.json")); aln = dict(CC.read_fasta(os.path.join(WORK, "aln", frame + ".fa")))
mods = jload(os.path.join(WORK, "aln", frame + "_modules.json"))["modules"] if os.path.exists(os.path.join(WORK, "aln", frame + "_modules.json")) else jload(os.path.join(WORK, "modules_%s.json" % frame))["modules"]
core = [(a, b) for n, a, b in mods if n in ("DI", "DII", "DIII")]
cols = np.concatenate([np.arange(a - 1, b) for a, b in core])
M = np.array([[ord(aln[r["id"]][c]) for c in cols] for r in rows], dtype=np.int16); gap = ord("-")
valid = M != gap; n = len(rows); D = np.zeros((n, n), dtype=np.float32)
for i in range(n):
    both = valid & valid[i]; same = (M == M[i]) & both; nb = both.sum(1)
    D[i] = np.where(nb >= 100, 1.0 - same.sum(1) / np.maximum(nb, 1), 1.0)
fam = [r["family"] for r in rows]; fams = sorted(set(fam), key=lambda f: int(f[3:]) if f[3:].isdigit() else 999); fi = {f: [k for k in range(n) if fam[k] == f] for f in fams}
FD = np.zeros((len(fams), len(fams)))
for a, fa in enumerate(fams):
    for b, fb in enumerate(fams): FD[a, b] = D[np.ix_(fi[fa], fi[fb])].mean() if a != b else D[np.ix_(fi[fa], fi[fa])].mean() * 0.0
# top-down: repeatedly split the cluster with most rows (2-way cut of its own average-linkage tree) when the smaller part has >= MIN_PART rows, until N_CLASS - 1 clusters; the families left over form 'other'
MIN_PART = 10; fidx = {f: k for k, f in enumerate(fams)}; rowsof = lambda fl: sum(len(fi[f]) for f in fl)


def split(cl):
    """cut the average-linkage tree of the families of cl into k = 2..12 parts; small parts (< MIN_PART rows) are merged into the nearest part; first k that leaves >= 2 parts"""
    ix = [fidx[f] for f in cl]; sub = FD[np.ix_(ix, ix)]; Zs = linkage(squareform(sub, checks=False), "average")
    for k in range(2, min(12, len(cl)) + 1):
        lab2 = fcluster(Zs, k, criterion="maxclust"); parts = {v: [f for f, l in zip(cl, lab2) if l == v] for v in set(lab2)}
        big = {v: p for v, p in parts.items() if rowsof(p) >= MIN_PART}
        if len(big) >= 2:
            for v, p in parts.items():
                if v in big: continue
                near = min(big, key=lambda w: np.mean([FD[fidx[a_], fidx[b_]] for a_ in p for b_ in big[w]])); big[near] = big[near] + p
            return list(big.values())
    return None


clusters = [list(fams)]
while len(clusters) < N_CLASS:
    for cl in sorted(clusters, key=lambda c: -rowsof(c)):
        parts = split(cl) if len(cl) >= 2 else None
        if parts and len(clusters) - 1 + len(parts) <= N_CLASS + 2:
            clusters.remove(cl); clusters += parts; break
    else:
        break
# the dominant lineage (Cry1, Cry3, Cry7, Cry8, Cry9 ... equally distant from each other, mean core p-distance 0.62-0.67) cannot be split by the tree of family distances: split it by primary rank
# (Cry1 | Cry3, Cry7, Cry8 | the rest) and pool the clusters of < 25 rows into 'other small families'
big = next(c for c in clusters if "Cry1" in c); clusters.remove(big)
clusters += [["Cry1"], [f for f in big if f in ("Cry3", "Cry7", "Cry8")], [f for f in big if f not in ("Cry1", "Cry3", "Cry7", "Cry8")]]
small = [c for c in clusters if rowsof(c) < 25]; clusters = [c for c in clusters if rowsof(c) >= 25] + ([sum(small, [])] if small else [])
groups = {k: cl for k, cl in enumerate(clusters)}
out = {}; names = {}
for l, fl in sorted(groups.items(), key=lambda kv: -sum(len(fi[f]) for f in kv[1])):
    top = sorted(fl, key=lambda f: -len(fi[f]))[:3]; nm = " / ".join(top) + (" +" + str(len(fl) - 3) if len(fl) > 3 else "")
    if len(fl) == 1: nm = fl[0]
    if small and sorted(fl) == sorted(sum(small, [])): nm = "other small families"
    names[l] = nm
    for f in fl: out[f] = nm
    print("%-26s families %2d rows %3d : %s" % (nm, len(fl), sum(len(fi[f]) for f in fl), " ".join(sorted(fl, key=lambda f: int(f[3:]) if f[3:].isdigit() else 999))))
if not os.environ.get("CRY_DRY"): jdump(out, os.path.join(WORK, "classes.json"), indent=0)
