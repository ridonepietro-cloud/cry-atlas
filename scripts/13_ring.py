"""RING 4.0 on every structure file of the Cry rows (AFDB PDB models + experimental mmCIF), 7 in parallel in WSL with the local RING copy; resumable.
Input  RAW/afdb/*.pdb, RAW/pdb/*.cif ; output RAW/ring_out/<file>_ringEdges (edge files stay outside the project folder: ~1 MB per model; the contact records derived from them are in the pack).
Run: py -I 13_ring.py"""
import os, sys, subprocess, glob, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
OUT = os.path.join(RAW, "ring_out"); os.makedirs(OUT, exist_ok=True)
files = sorted(glob.glob(os.path.join(RAW, "afdb", "*.pdb")) + glob.glob(os.path.join(RAW, "pdb", "*.cif")))
todo = [f for f in files if not (os.path.exists(os.path.join(OUT, os.path.basename(f) + "_ringEdges")) and os.path.getsize(os.path.join(OUT, os.path.basename(f) + "_ringEdges")) > 50)]
print("structures", len(files), "to run", len(todo), flush=True)
sh = os.path.join(WORK, "run_ring.sh"); lst = os.path.join(WORK, "ring_todo.txt")
open(lst, "w", newline="\n").write("\n".join(wsl(f) for f in todo) + "\n")
open(sh, "w", newline="\n").write('''#!/bin/bash
# RING 4.0 (local copy), one structure per call, 7 in parallel; the edge file is written next to nothing else: stage in $H, copy back
set -u
RING="$HOME/tools/ring_local/out/bin/ring"; H="$HOME/cry_ring"; OUTD="%s"
mkdir -p "$H/in" "$H/out"
run_one() {
  f="$1"; b=$(basename "$f"); RING="$2"; H="$3"; OUTD="$4"
  cp "$f" "$H/in/$b"; "$RING" -i "$H/in/$b" --out_dir "$H/out" > /dev/null 2>&1
  if [ -s "$H/out/${b}_ringEdges" ]; then cp "$H/out/${b}_ringEdges" "$OUTD/"; fi
  rm -f "$H/in/$b" "$H/out/${b}_ringNodes" "$H/out/${b}_ringEdges"
}
export -f run_one
cat "%s" | xargs -P 7 -I{} bash -c 'run_one "$@"' _ {} "$RING" "$H" "$OUTD"
echo done
''' % (wsl(OUT), wsl(lst)))
t0 = time.time()
subprocess.run(["wsl", "-e", "bash", wsl(sh)], check=True)
got = sum(1 for f in files if os.path.exists(os.path.join(OUT, os.path.basename(f) + "_ringEdges")))
print("edge files %d of %d structures (%.0f s)" % (got, len(files), time.time() - t0))
print("missing:", [os.path.basename(f) for f in files if not os.path.exists(os.path.join(OUT, os.path.basename(f) + "_ringEdges"))][:20])
