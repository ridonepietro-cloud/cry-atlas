"""Name-order trees of the Cry pack -> site/data/cry/trees.json ["name_order"]  (one tree per alignment frame: structure-based and sequence-based).

What: a tree CALCULATED FROM THE ENTRIES IN NAME ORDER, offered in the tree window next to the ML tree. Its tips read top to bottom in the Name order of the page, so with Order = Name tips and rows are level
(no connecting lines). Such a tree can only be an ORDER-CONSTRAINED tree: average-linkage clustering (UPGMA) in which only clusters that are neighbours in the Name list may be joined.
It is a distance tree, not a maximum-likelihood tree: use the ML trees for relationships, this one to see how far the names follow the sequences / structures.

Distances: for every pair of entries the fraction p of differing residues over the columns where both have a residue (columns occupied by >= 20 % of the entries, as for the ML trees of
the full toxin; pairs sharing < 50 such columns get the largest distance), corrected as Kimura (1983) protein distance  d = -ln(1 - p - 0.2 p^2)  (p capped at 0.85).
Clustering: start with one cluster per entry in Name order; repeatedly join the two ADJACENT clusters with the smallest mean pairwise distance; node height = max(mean / 2, heights of its two children)
(the constraint can make a later join smaller than an earlier one; the height is then raised to keep the tree drawable: the number of such joins is stored).
Name order = the order of the page (label, case-insensitive, character by character; the page uses the order stored here, so tips and rows cannot disagree between browser locales).
Run: py 28_name_order_trees.py     (NOT -I if the unconstrained reference correlation should be stored: it needs scipy of the user site)
Reads site/data/cry/{meta,rows,aln_*,trees}.json and adds / replaces trees.json["name_order"]. 20_trees.py rewrites trees.json from scratch: run this script after it."""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import numpy as np
D = os.path.join(SITE, "data", "cry"); MIN_OCC = 0.20; MIN_SHARED = 50; P_CAP = 0.85
meta = jload(os.path.join(D, "meta.json")); rows = jload(os.path.join(D, "rows.json")); trees = jload(os.path.join(D, "trees.json"))
order = [r["id"] for r in sorted(rows, key=lambda r: r["label"].casefold())]
AA = "ACDEFGHIKLMNPQRSTVWY"


def distances(seqs):
    """Kimura protein distances between all entries (in `order`) over the columns occupied by >= MIN_OCC of the entries"""
    n = len(order)
    M = np.array([list(seqs[i]) for i in order], dtype="U1")
    occ = (M != "-") & (M != "X") & (M != "?")
    keep = occ.sum(0) >= max(1, MIN_OCC * n); M = M[:, keep]; occ = occ[:, keep]; ncol = int(keep.sum())
    O = occ.astype(np.float32); shared = O @ O.T
    ident = np.zeros((n, n), dtype=np.float32)
    for a in AA:
        X = (M == a).astype(np.float32); ident += X @ X.T
    with np.errstate(divide="ignore", invalid="ignore"):
        p = 1.0 - ident / shared
    bad = shared < MIN_SHARED; p = np.where(bad, P_CAP, np.minimum(p, P_CAP)); p = np.nan_to_num(p, nan=P_CAP)
    d = -np.log(1.0 - p - 0.2 * p * p); np.fill_diagonal(d, 0.0)
    return d.astype(np.float64), ncol, int(bad.sum() // 2)


def constrained_upgma(d):
    """order-constrained average linkage -> (newick, cophenetic matrix, number of height corrections)"""
    n = d.shape[0]; coph = np.zeros_like(d)
    cl = [{"lo": i, "hi": i, "h": 0.0, "nw": order[i]} for i in range(n)]       # clusters in list order
    def mean(a, b): return float(d[a["lo"]:a["hi"] + 1, b["lo"]:b["hi"] + 1].mean())
    link = [mean(cl[i], cl[i + 1]) for i in range(n - 1)]; fixes = 0
    while len(cl) > 1:
        k = int(np.argmin(link)); a, b = cl[k], cl[k + 1]; m = link[k]; h = m / 2.0
        if h < max(a["h"], b["h"]): h = max(a["h"], b["h"]); fixes += 1
        coph[a["lo"]:a["hi"] + 1, b["lo"]:b["hi"] + 1] = 2 * h; coph[b["lo"]:b["hi"] + 1, a["lo"]:a["hi"] + 1] = 2 * h
        nw = "(%s:%.6f,%s:%.6f)" % (a["nw"], h - a["h"], b["nw"], h - b["h"])
        c = {"lo": a["lo"], "hi": b["hi"], "h": h, "nw": nw}
        cl[k:k + 2] = [c]; del link[k]                    # the join of a and b is gone; link[k - 1] and link[k] now touch the new cluster
        if k > 0: link[k - 1] = mean(cl[k - 1], c)
        if k < len(cl) - 1: link[k] = mean(c, cl[k + 1])
    return cl[0]["nw"] + ";", coph, fixes


def corr(a, b):
    iu = np.triu_indices(a.shape[0], 1)
    return float(np.corrcoef(a[iu], b[iu])[0, 1])


out = {"order": order, "method": "UPGMA (average linkage) constrained to the Name order: only clusters that are neighbours in the Name list are joined; Kimura protein distances over the alignment columns occupied by >= 20 % of the entries"}
for fr in [f["id"] for f in meta["frames"]]:
    seqs = jload(os.path.join(D, next(f["file"] for f in meta["frames"] if f["id"] == fr)))["seqs"]
    d, ncol, nbad = distances(seqs)
    nw, coph, fixes = constrained_upgma(d)
    info = {"newick": nw, "n_taxa": len(order), "n_columns": ncol, "pairs_few_shared_columns": nbad, "height_corrections": fixes, "cophenetic_r": round(corr(coph, d), 4),
            "label": "Cry toxin · name-order tree · %s" % ("structure-based" if fr == "struct" else "sequence-based")}
    try:                                                   # reference: the same clustering without the order constraint
        from scipy.cluster.hierarchy import linkage, cophenet
        from scipy.spatial.distance import squareform
        Z = linkage(squareform(d, checks=False), "average"); info["cophenetic_r_unconstrained"] = round(float(np.corrcoef(cophenet(Z), squareform(d, checks=False))[0, 1]), 4)
    except Exception:
        info["cophenetic_r_unconstrained"] = None
    out[fr] = info
    print(fr, {k: v for k, v in info.items() if k != "newick"}, flush=True)
trees["name_order"] = out
jdump(trees, os.path.join(D, "trees.json"), separators=(",", ":"))
print("trees.json updated:", sorted(trees.keys()), os.path.getsize(os.path.join(D, "trees.json")) // 1024, "KB")
