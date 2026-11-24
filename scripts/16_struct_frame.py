"""Structure-based frame of the Cry Atlas: FoldMason `easy-msa --report-mode 1` of the AlphaFold DB models whose residues equal the row sequence (exact-sequence models),
every other row (no exact model, donor-model rows, rows without a model) is added with `mafft --add` to the amino-acid alignment of FoldMason (same mechanism as the other packs).
A frame row must degap to the row sequence (checked).  Output work/aln/struct.fa, work/aln/struct_report.json.   Run: py -I 16_struct_frame.py   (needs work/rows.json from 14_assemble_rows.py)"""
import os, sys, json, time, shutil, subprocess, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
AL = os.path.join(WORK, "aln"); os.makedirs(AL, exist_ok=True)
THREADS = 6
MAFFT_HEAD = "export MAFFT_BINARIES=$HOME/tools/mafft-linux64/mafftdir/libexec\nM=$HOME/tools/mafft-linux64/mafftdir/bin/mafft\n"
rows = jload(os.path.join(WORK, "rows.json")); seq = {r["id"]: r["seq"] for r in rows}
t0 = time.time(); pdbs = os.path.join(AL, "fm_pdb"); os.makedirs(pdbs, exist_ok=True)
for fn in os.listdir(pdbs): os.remove(os.path.join(pdbs, fn))
cand = []
for r in rows:
    m = next((m for m in r["models"] if m["kind"] == "AF2" and "donor" not in m and os.path.exists(m["struct"])), None)
    if m:
        shutil.copyfile(m["struct"], os.path.join(pdbs, r["id"] + ".pdb")); cand.append(r["id"])
print("FoldMason input models:", len(cand), flush=True)
sh = ["#!/bin/bash", 'FM="$HOME/tools/foldmason/bin/foldmason"', 'SRC="%s"' % wsl(AL), 'W="$HOME/fm_cry"', 'rm -rf "$W"; mkdir -p "$W/out"; cp -r "$SRC/fm_pdb" "$W/pdb"; cd "$W"',
      'n=$(ls pdb | wc -l); t0=$(date +%s); "$FM" easy-msa pdb out/core tmp --report-mode 1 --threads 6 > fm.log 2>&1; rc=$?',
      'echo "foldmason rc=$rc models=$n seconds=$(( $(date +%s) - t0 ))" > "$SRC/fm_done.txt"; grep -m1 "Average MSA LDDT" fm.log >> "$SRC/fm_done.txt"',
      'cp out/core_aa.fa "$SRC/fm_aa.fa"; cp out/core_3di.fa "$SRC/fm_3di.fa"; cp fm.log "$SRC/fm.log"; cp out/core.nw "$SRC/fm.nw" 2>/dev/null; rm -rf tmp']
p = os.path.join(AL, "run_fm.sh"); open(p, "w", newline="\n").write("\n".join(sh) + "\n")
subprocess.run(["wsl", "-e", "bash", wsl(p)], check=True)
print(open(os.path.join(AL, "fm_done.txt")).read(), flush=True)
fm = collections.OrderedDict((h.split()[0].replace(".pdb", ""), s.upper()) for h, s in CC.read_fasta(os.path.join(AL, "fm_aa.fa")))
ok = {i: s for i, s in fm.items() if i in seq and s.replace("-", "") == seq[i]}
bad = [i for i in fm if i not in ok]
print("FoldMason: %d models, %d rows equal their sequence, %d differ" % (len(fm), len(ok), len(bad)), flush=True)
rest = [(r["id"], r["seq"]) for r in rows if r["id"] not in ok]
with open(os.path.join(AL, "fm_ok.fa"), "w", newline="\n") as f:
    for i, s in ok.items(): f.write(">%s\n%s\n" % (i, s))
with open(os.path.join(AL, "struct_rest.faa"), "w", newline="\n") as f:
    for i, s in rest: f.write(">%s\n%s\n" % (i, s))
sh = [MAFFT_HEAD, 'cd "%s"' % wsl(AL), "$M --add struct_rest.faa --thread %d fm_ok.fa > struct_all.aln 2> struct_add.log" % THREADS]
p = os.path.join(AL, "run_struct_add.sh"); open(p, "w", newline="\n").write("\n".join(sh) + "\n")
subprocess.run(["wsl", "-e", "bash", wsl(p)], check=True)
aln = collections.OrderedDict((h.split()[0], s.upper()) for h, s in CC.read_fasta(os.path.join(AL, "struct_all.aln")))
for i, s in seq.items(): assert aln[i].replace("-", "") == s, "struct frame row %s differs" % i
with open(os.path.join(AL, "struct.fa"), "w", newline="\n") as f:
    for i, s in aln.items(): f.write(">%s\n%s\n" % (i, s))
jdump({"models": len(fm), "foldmason_rows": len(ok), "added_with_mafft": len(rest), "width": len(next(iter(aln.values()))), "seconds": round(time.time() - t0)}, os.path.join(AL, "struct_report.json"), indent=1)
print("struct frame width", len(next(iter(aln.values()))), "seconds", round(time.time() - t0))
