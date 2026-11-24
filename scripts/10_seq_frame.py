"""Sequence-based frame of the Cry Atlas (all 670 rows, full-length sequences).
MAFFT L-INS-i (--localpair --maxiterate 100) of cluster representatives (greedy clustering by shared 5-mers at CL identity), then `mafft --add` of the other rows (as the FliC pack).
Output work/aln/seq.fa (degaps to every row sequence; checked), work/aln/seq_order.json (MAFFT --reorder guide-tree order), work/aln/seq_report.json.
Run: py -I 10_seq_frame.py [--dry]"""
import os, sys, re, json, time, collections, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
AL = os.path.join(WORK, "aln"); os.makedirs(AL, exist_ok=True)
CL = float(os.environ.get("CRY_CL", "0.75")); THREADS = 6
MAFFT_HEAD = "export MAFFT_BINARIES=$HOME/tools/mafft-linux64/mafftdir/libexec\nM=$HOME/tools/mafft-linux64/mafftdir/bin/mafft\n"


def kset(s, k=5): return {s[i:i + k] for i in range(len(s) - k + 1)}


def cluster(items, ident):
    thr = 0.8 * ident ** 5
    reps, ks, of = [], {}, {}
    for i, s in sorted(items, key=lambda t: -len(t[1])):
        k = kset(s); hit = None
        for r in reps:
            if len(k & ks[r]) / max(1, min(len(k), len(ks[r]))) >= thr: hit = r; break
        if hit is None: reps.append(i); ks[i] = k; of[i] = i
        else: of[i] = hit
    return reps, of


rows = jload(os.path.join(WORK, "rows_seed.json")); items = [(r["id"], r["seq"]) for r in rows]
reps, of = cluster(items, CL)
print("rows", len(items), "clusters at", CL, "->", len(reps), "representatives; cluster sizes", collections.Counter(collections.Counter(of.values()).values()).most_common(6), flush=True)
if "--dry" in sys.argv: sys.exit()
seq = dict(items); t0 = time.time()
with open(os.path.join(AL, "seq_reps.faa"), "w", newline="\n") as f:
    for i in reps: f.write(">%s\n%s\n" % (i, seq[i]))
with open(os.path.join(AL, "seq_rest.faa"), "w", newline="\n") as f:
    for i, s in items:
        if i not in set(reps): f.write(">%s\n%s\n" % (i, s))
sh = [MAFFT_HEAD, 'cd "%s"' % wsl(AL), "$M --localpair --maxiterate 100 --reorder --thread %d seq_reps.faa > seq_reps.aln 2> seq_reps.log" % THREADS,
      "$M --add seq_rest.faa --reorder --thread %d seq_reps.aln > seq_all.aln 2> seq_add.log" % THREADS]
p = os.path.join(AL, "run_seq.sh"); open(p, "w", newline="\n").write("\n".join(sh) + "\n")
subprocess.run(["wsl", "-e", "bash", wsl(p)], check=True)
aln = collections.OrderedDict((h.split()[0], s.upper()) for h, s in read_fasta(os.path.join(AL, "seq_all.aln")))
for i, s in items: assert aln[i].replace("-", "") == s, "frame row %s differs" % i
with open(os.path.join(AL, "seq.fa"), "w", newline="\n") as f:
    for i, s in aln.items(): f.write(">%s\n%s\n" % (i, s))
jdump({"order": list(aln.keys())}, os.path.join(AL, "seq_order.json"))
jdump({"reps": len(reps), "rest": len(items) - len(reps), "cl": CL, "width": len(next(iter(aln.values()))), "seconds": round(time.time() - t0)}, os.path.join(AL, "seq_report.json"), indent=1)
print("done: width", len(next(iter(aln.values()))), "seconds", round(time.time() - t0))
