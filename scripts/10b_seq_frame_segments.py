"""Sequence-based frame of the Cry Atlas, SEGMENT-WISE (replaces the single-pass frame of 10_seq_frame.py, which stretched Domain I over 610 columns because local alignment staggers the non-homologous
N- and C-terminal regions).
 1. segments of every row = slices of its sequence by the module of the column of each residue in the first-pass frame (work/aln/seq.fa + work/modules_seq.json: N-term | DI | DII | DIII | C-term; a local
    alignment keeps the residue order, so the slices are contiguous);
 2. every segment is aligned on its own: MAFFT L-INS-i (--localpair --maxiterate 100) of cluster representatives (>= MINLEN residues; greedy 5-mer clustering at CL), the other rows added with `mafft --add`;
 3. the segment alignments are joined in order; a row without residues in a segment has gaps there.   (the Atlas joins the A and B subunits the same way)
Output work/aln/seqseg.fa (degaps to every row sequence; checked), work/aln/seqseg_modules.json (module column ranges = segment widths), work/aln/seqseg_order.json, report.   Run: py -I 10b_seq_frame_segments.py [--dry]"""
import os, sys, re, json, time, collections, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
AL = os.path.join(WORK, "aln", "seg"); os.makedirs(AL, exist_ok=True)
CL = float(os.environ.get("CRY_CL", "0.85")); MINLEN = 30; THREADS = 6
MAFFT_HEAD = "export MAFFT_BINARIES=$HOME/tools/mafft-linux64/mafftdir/libexec\nM=$HOME/tools/mafft-linux64/mafftdir/bin/mafft\n"
rows = jload(os.path.join(WORK, "rows_seed.json")); first = dict(CC.read_fasta(os.path.join(WORK, "aln", "seq.fa"))); mods = jload(os.path.join(WORK, "modules_seq.json"))["modules"]
NAMES = [m[0] for m in mods]; col2mod = {}
for name, a, b in mods:
    for c in range(a, b + 1): col2mod[c] = name
seg = {n: {} for n in NAMES}
for r in rows:
    cols = CC.col_index(first[r["id"]]); cur = {n: [] for n in NAMES}
    for i, c in enumerate(cols): cur[col2mod[c]].append(i)
    last = -1
    for n in NAMES:
        ix = cur[n]
        if not ix: continue
        assert ix == list(range(ix[0], ix[-1] + 1)) and ix[0] > last, "row %s: segment %s is not a contiguous slice" % (r["id"], n)
        last = ix[-1]; seg[n][r["id"]] = r["seq"][ix[0]:ix[-1] + 1]
for n in NAMES:
    L = [len(s) for s in seg[n].values()]
    print("%-7s rows %3d  length median %d max %d  (rows < %d aa: %d)" % (n, len(L), sorted(L)[len(L) // 2], max(L), MINLEN, sum(1 for x in L if x < MINLEN)), flush=True)


def kset(s, k=5): return {s[i:i + k] for i in range(len(s) - k + 1)}


def cluster(items, ident):
    thr = 0.8 * ident ** 5; reps, ks = [], {}
    for i, s in sorted(items, key=lambda t: -len(t[1])):
        k = kset(s)
        if not any(len(k & ks[r]) / max(1, min(len(k), len(ks[r]))) >= thr for r in reps): reps.append(i); ks[i] = k
    return reps


plan = {}
for n in NAMES:
    items = [(i, s) for i, s in seg[n].items()]; big = [(i, s) for i, s in items if len(s) >= MINLEN]
    reps = cluster(big, CL) if len(big) > 3 else [i for i, _ in big]
    plan[n] = (reps, [i for i, _ in items if i not in set(reps)]); print(n, "representatives", len(reps), "added", len(plan[n][1]), flush=True)
if "--dry" in sys.argv: sys.exit()
t0 = time.time(); sh = [MAFFT_HEAD, 'cd "%s"' % wsl(AL)]
for k, n in enumerate(NAMES):
    tag = "s%d" % k; reps, rest = plan[n]
    with open(os.path.join(AL, tag + "_reps.faa"), "w", newline="\n") as f:
        for i in reps: f.write(">%s\n%s\n" % (i, seg[n][i]))
    with open(os.path.join(AL, tag + "_rest.faa"), "w", newline="\n") as f:
        for i in rest: f.write(">%s\n%s\n" % (i, seg[n][i]))
    if len(reps) >= 2:
        sh.append("$M --localpair --maxiterate 100 --reorder --thread %d %s_reps.faa > %s_reps.aln 2> %s_reps.log" % (THREADS, tag, tag, tag))
    else:
        sh.append("cp %s_reps.faa %s_reps.aln" % (tag, tag))
    if rest:
        sh.append("$M --add %s_rest.faa --reorder --thread %d %s_reps.aln > %s_all.aln 2> %s_add.log" % (tag, THREADS, tag, tag, tag))
    else:
        sh.append("cp %s_reps.aln %s_all.aln" % (tag, tag))
p = os.path.join(AL, "run_seg.sh"); open(p, "w", newline="\n").write("\n".join(sh) + "\n")
subprocess.run(["wsl", "-e", "bash", wsl(p)], check=True)
aligned = {}; widths = {}
for k, n in enumerate(NAMES):
    a = collections.OrderedDict((h.split()[0], s.upper()) for h, s in CC.read_fasta(os.path.join(AL, "s%d_all.aln" % k)))
    for i, s in seg[n].items(): assert a[i].replace("-", "") == s, "segment %s row %s differs" % (n, i)
    aligned[n] = a; widths[n] = len(next(iter(a.values()))) if a else 0
frame = collections.OrderedDict(); order = []
for r in rows:
    frame[r["id"]] = "".join(aligned[n].get(r["id"], "-" * widths[n]) for n in NAMES)
    assert frame[r["id"]].replace("-", "") == r["seq"], "joined row %s differs" % r["id"]
modules = []; c = 1
for n in NAMES:
    if widths[n]: modules.append([n, c, c + widths[n] - 1]); c += widths[n]
with open(os.path.join(WORK, "aln", "seqseg.fa"), "w", newline="\n") as f:
    for i, s in frame.items(): f.write(">%s\n%s\n" % (i, s))
# guide-tree order: order of the DII segment alignment (largest homologous domain), then rows without DII in id order
ordr = list(aligned["DII"].keys()); ordr += [r["id"] for r in rows if r["id"] not in set(ordr)]
jdump({"order": ordr}, os.path.join(WORK, "aln", "seqseg_order.json")); jdump({"modules": modules, "widths": widths}, os.path.join(WORK, "aln", "seqseg_modules.json"), indent=1)
jdump({"cl": CL, "reps": {n: len(plan[n][0]) for n in NAMES}, "width": c - 1, "seconds": round(time.time() - t0)}, os.path.join(WORK, "aln", "seqseg_report.json"), indent=1)
print("frame width", c - 1, "modules", modules, "seconds", round(time.time() - t0))
