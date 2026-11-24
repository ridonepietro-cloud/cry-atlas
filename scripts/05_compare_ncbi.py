"""Compare BPPRC sequences with the NCBI GenPept sequence of the same accession -> work/ncbi/compare.tsv and a printed summary."""
import os, sys, re, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
tab = list(csv.DictReader(open(os.path.join(DATA, "bpprc_cry_table.tsv"), encoding="utf8"), delimiter="\t"))
seqs = jload(os.path.join(DATA, "cry_seq_clean.json"))["seqs"]; gp = jload(os.path.join(WORK, "ncbi", "genpept.json"))
from Bio.Align import PairwiseAligner
al = PairwiseAligner(); al.mode = "global"; al.match_score = 1; al.mismatch_score = 0; al.open_gap_score = -1; al.extend_gap_score = -0.5; al.end_gap_score = 0
cnt = collections.Counter(); rows = []
for r in tab:
    n, a = r["name"], r["accession"]; s = seqs.get(n); g = gp.get(a)
    if s is None: cnt["no BPPRC sequence"] += 1; continue
    if not g or g.get("missing"): cnt["no NCBI protein record"] += 1; rows.append((n, a, len(s), "", "no NCBI protein record", "")); continue
    t = g["seq"]
    if t == s: cnt["identical"] += 1; rows.append((n, a, len(s), len(t), "identical", "")); continue
    if t.replace("X", "") == s.replace("X", ""): cnt["identical except X"] += 1
    aln = al.align(s, t)[0]
    same = sum(1 for (a0, a1), (b0, b1) in zip(*aln.aligned) for k in range(a1 - a0) if s[a0 + k] == t[b0 + k])
    cov = sum(a1 - a0 for (a0, a1), _ in zip(*aln.aligned))
    kind = "different length" if len(s) != len(t) else "same length, substitutions"
    cnt[kind] += 1
    rows.append((n, a, len(s), len(t), kind, "aligned %d, identical %d (%.1f%%) ; ncbi has %d, bpprc %d" % (cov, same, 100.0 * same / max(1, cov), len(t), len(s))))
print(cnt)
with open(os.path.join(WORK, "ncbi", "compare.tsv"), "w", encoding="utf8", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["name", "accession", "len_bpprc", "len_ncbi", "result", "detail"]); w.writerows(rows)
bad = [x for x in rows if x[4] not in ("identical", "no NCBI protein record")]
print(len(bad), "differing"); 
for x in bad[:60]: print(x)
