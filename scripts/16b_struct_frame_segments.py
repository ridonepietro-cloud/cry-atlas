"""Structure-based frame of the Cry Atlas, SEGMENT-WISE (replaces the whole-model FoldMason frame of 16_struct_frame.py: average LDDT 0.26 and 5,686 columns, because the protoxin halves of the
full-length models differ).  Per module (N-term | DI | DII | DIII | C-term; slices of the sequence from the segment frame work/aln/seqseg.fa):
 1. the AlphaFold DB models (exact-sequence models) are cut to the residues of the segment (>= MINLEN residues) -> FoldMason `easy-msa --report-mode 1` (amino-acid alignment of the structures);
 2. the other rows of the segment (no exact model, donor-model rows, rows without a model, short segments) are added with `mafft --add`;
 3. the segment alignments are joined in order.   A frame row must degap to the row sequence (checked).
Output work/aln/structseg.fa, work/aln/structseg_modules.json, work/aln/structseg_report.json.   Run: py -I 16b_struct_frame_segments.py"""
import os, sys, json, time, shutil, subprocess, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
AL = os.path.join(WORK, "aln", "structseg"); os.makedirs(AL, exist_ok=True)
THREADS = 6; MINLEN = 30
MAFFT_HEAD = "export MAFFT_BINARIES=$HOME/tools/mafft-linux64/mafftdir/libexec\nM=$HOME/tools/mafft-linux64/mafftdir/bin/mafft\n"
rows = jload(os.path.join(WORK, "rows.json")); seq = {r["id"]: r["seq"] for r in rows}
segf = dict(CC.read_fasta(os.path.join(WORK, "aln", "seqseg.fa"))); segm = jload(os.path.join(WORK, "aln", "seqseg_modules.json"))["modules"]; NAMES = [m[0] for m in segm]
col2mod = {}
for n, a, b in segm:
    for c in range(a, b + 1): col2mod[c] = n
rng = {n: {} for n in NAMES}                                           # row -> (first, last) 0-based residue index of the segment
for r in rows:
    cur = collections.defaultdict(list)
    for i, c in enumerate(CC.col_index(segf[r["id"]])): cur[col2mod[c]].append(i)
    for n, ix in cur.items(): rng[n][r["id"]] = (ix[0], ix[-1])
model = {}
for r in rows:
    m = next((m for m in r["models"] if m["kind"] == "AF2" and "donor" not in m and os.path.exists(m["struct"])), None)
    if m: model[r["id"]] = m["struct"]
t0 = time.time(); sh = ["#!/bin/bash", 'FM="$HOME/tools/foldmason/bin/foldmason"', 'SRC="%s"' % wsl(AL), "set -u"]; plan = {}
for k, n in enumerate(NAMES):
    d = os.path.join(AL, "pdb_%d" % k); os.makedirs(d, exist_ok=True)
    for fn in os.listdir(d): os.remove(os.path.join(d, fn))
    cand = []
    for rid, (a, b) in rng[n].items():
        if rid in model and b - a + 1 >= MINLEN:
            keep = []
            for ln in open(model[rid], encoding="utf-8", errors="replace"):
                if ln.startswith("ATOM") and a + 1 <= int(ln[22:26]) <= b + 1: keep.append(ln)
            open(os.path.join(d, rid + ".pdb"), "w", newline="\n").write("".join(keep) + "END\n"); cand.append(rid)
    plan[n] = cand; print(n, "segment structures:", len(cand), "of", len(rng[n]), "rows", flush=True)
    sh.append('W="$HOME/fm_cry_seg/%d"; rm -rf "$W"; mkdir -p "$W/out"; cp -r "$SRC/pdb_%d" "$W/pdb"; cd "$W"; "$FM" easy-msa pdb out/core tmp --report-mode 1 --threads %d > fm.log 2>&1; echo "seg %d rc=$?" >> "$SRC/fm_done.txt"; grep -m1 "Average MSA LDDT" fm.log >> "$SRC/fm_done.txt"; cp out/core_aa.fa "$SRC/fm_aa_%d.fa"; rm -rf tmp' % (k, k, THREADS, k, k))
open(os.path.join(AL, "fm_done.txt"), "w").write("")
p = os.path.join(AL, "run_fm.sh"); open(p, "w", newline="\n").write("\n".join(sh) + "\n")
subprocess.run(["wsl", "-e", "bash", wsl(p)], check=True)
print(open(os.path.join(AL, "fm_done.txt")).read(), flush=True)
sh = [MAFFT_HEAD, 'cd "%s"' % wsl(AL)]; seqs = {}; ok = {}
for k, n in enumerate(NAMES):
    seqs[n] = {rid: seq[rid][a:b + 1] for rid, (a, b) in rng[n].items()}
    fm = collections.OrderedDict((h.split()[0].replace(".pdb", ""), s.upper()) for h, s in CC.read_fasta(os.path.join(AL, "fm_aa_%d.fa" % k)))
    ok[n] = {i: s for i, s in fm.items() if i in seqs[n] and s.replace("-", "") == seqs[n][i]}
    bad = [i for i in fm if i not in ok[n]]; rest = [(i, s) for i, s in seqs[n].items() if i not in ok[n]]
    print(n, "FoldMason rows equal to their segment:", len(ok[n]), "differ:", len(bad), "added with MAFFT:", len(rest), flush=True)
    with open(os.path.join(AL, "fm_ok_%d.fa" % k), "w", newline="\n") as f:
        for i, s in ok[n].items(): f.write(">%s\n%s\n" % (i, s))
    with open(os.path.join(AL, "rest_%d.faa" % k), "w", newline="\n") as f:
        for i, s in rest: f.write(">%s\n%s\n" % (i, s))
    sh.append("$M --add rest_%d.faa --thread %d fm_ok_%d.fa > all_%d.aln 2> add_%d.log" % (k, THREADS, k, k, k) if rest else "cp fm_ok_%d.fa all_%d.aln" % (k, k))
p = os.path.join(AL, "run_add.sh"); open(p, "w", newline="\n").write("\n".join(sh) + "\n")
subprocess.run(["wsl", "-e", "bash", wsl(p)], check=True)
aligned = {}; widths = {}
for k, n in enumerate(NAMES):
    a = collections.OrderedDict((h.split()[0], s.upper()) for h, s in CC.read_fasta(os.path.join(AL, "all_%d.aln" % k)))
    for i, s in seqs[n].items(): assert a[i].replace("-", "") == s, "segment %s row %s differs" % (n, i)
    aligned[n] = a; widths[n] = len(next(iter(a.values())))
frame = collections.OrderedDict()
for r in rows:
    frame[r["id"]] = "".join(aligned[n].get(r["id"], "-" * widths[n]) for n in NAMES); assert frame[r["id"]].replace("-", "") == r["seq"], "joined row %s differs" % r["id"]
modules = []; c = 1
for n in NAMES:
    modules.append([n, c, c + widths[n] - 1]); c += widths[n]
with open(os.path.join(WORK, "aln", "structseg.fa"), "w", newline="\n") as f:
    for i, s in frame.items(): f.write(">%s\n%s\n" % (i, s))
jdump({"modules": modules, "widths": widths}, os.path.join(WORK, "aln", "structseg_modules.json"), indent=1)
jdump({"models": {n: len(ok[n]) for n in NAMES}, "width": c - 1, "seconds": round(time.time() - t0), "lddt": open(os.path.join(AL, "fm_done.txt")).read()}, os.path.join(WORK, "aln", "structseg_report.json"), indent=1)
print("structure segment frame width", c - 1, "modules", modules, "seconds", round(time.time() - t0))
