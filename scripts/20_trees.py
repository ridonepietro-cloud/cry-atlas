"""Trees of the Cry pack -> site/data/cry/trees.json (generalisation of build_trees.py of the NODE-MSA Atlas for one subunit with several module trees).
Trees per frame: A = all columns occupied by >= MIN_OCC_A of the rows (the variable N- and C-terminal insertions are left out); DI, DII, DIII = the columns of that module, rows with >= MIN_RES residues of it
(domain trees for the tanglegram: Cry toxins are known for domain swaps).   IQ-TREE 2.3.6 (WSL), fixed LG+G4, --fast search, SH-aLRT support (1000 replicates), seed 20261124; results are cached by a hash.
Run: py -I 20_trees.py [frames=struct,seq] [sets=A,DI,DII,DIII] [threads=3] [parallel=2]"""
import os, sys, json, hashlib, subprocess, re, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
D = os.path.join(SITE, "data", "cry"); W = os.path.join(WORK, "cry", "trees"); os.makedirs(W, exist_ok=True)
args = dict(a.split("=") for a in sys.argv[1:] if "=" in a)
FRAMES = args.get("frames", "struct,seq").split(","); SETS = args.get("sets", "A,DI,DII,DIII").split(","); THREADS = int(args.get("threads", 3)); PAR = int(args.get("parallel", 2))
IQ = "$HOME/tools/iqtree-2.3.6-Linux-intel/bin/iqtree2"; SEED = 20261124; CMD = "-m LG+G4 --fast -alrt 1000"
MIN_OCC = {"A": 0.20, "DI": 0.10, "DII": 0.10, "DIII": 0.10}; MIN_RES = {"DI": 100, "DII": 100, "DIII": 70}
meta = jload(os.path.join(D, "meta.json")); mods = {m["name"]: m["cols"] for m in next(iter(meta["modules"].values()))}
LABEL = {"A": "full toxin", "DI": "Domain I", "DII": "Domain II", "DIII": "Domain III"}


def sets_of(frame):
    seqs = jload(os.path.join(D, next(f["file"] for f in meta["frames"] if f["id"] == frame)))["seqs"]; out = {}
    # module column ranges of THIS frame: the page stores them for the first frame only; recompute from the modules of the row residues
    res = jload(os.path.join(D, "residues.json"))
    modcols = {}
    for rid, g in seqs.items():
        mod = res[rid]["A"]["mod"]; import base64; mb = base64.b64decode(mod); k = 0
        for c, ch in enumerate(g):
            if ch != "-":
                modcols.setdefault(mb[k], set()).add(c); k += 1
    name_of = meta["module_ids"]                                 # module id (string) -> name
    cols_of = {name_of[str(m)]: sorted(c) for m, c in modcols.items() if str(m) in name_of}
    for s in SETS:
        cols = list(range(len(next(iter(seqs.values()))))) if s == "A" else cols_of.get(s, [])
        part = {rid: "".join(g[c] for c in cols) for rid, g in seqs.items()}
        if s != "A": part = {rid: x for rid, x in part.items() if len(x.replace("-", "")) >= MIN_RES[s]}
        else: part = {rid: x for rid, x in part.items() if x.replace("-", "")}
        need = max(1, MIN_OCC[s] * len(part)); keep = [c for c in range(len(next(iter(part.values())))) if sum(1 for x in part.values() if x[c] != "-") >= need]
        out[s] = {rid: "".join(x[c] for c in keep) for rid, x in part.items()}
    return out


jobs = {}
for fr in FRAMES:
    for s, aln in sets_of(fr).items():
        key = "%s_%s" % (fr, s); fa = os.path.join(W, key + ".fa")
        with open(fa, "w", newline="\n") as f:
            for rid in sorted(aln): f.write(">%s\n%s\n" % (rid, aln[rid]))
        sig = hashlib.sha256(open(fa, "rb").read() + (CMD + str(SEED)).encode()).hexdigest()
        tf = os.path.join(W, key + ".treefile"); stamp = os.path.join(W, key + ".sha")
        cached = os.path.exists(tf) and os.path.exists(stamp) and open(stamp).read().strip() == sig
        if not cached:
            for ext in (".treefile", ".iqtree", ".log", ".contree", ".mldist", ".splits.nex", ".uniqueseq.phy", ".model.gz", ".bionj", ".ckp.gz"):
                if os.path.exists(os.path.join(W, key + ext)): os.remove(os.path.join(W, key + ext))
            sh = os.path.join(W, "run_%s.sh" % key)
            open(sh, "wb").write(("#!/bin/bash\ncd \"%s\"\n%s -s %s.fa %s -T %d -seed %d --prefix %s -redo --quiet\n" % (wsl(W), IQ, key, CMD, THREADS, SEED, key)).encode())
        jobs[key] = {"aln": aln, "tf": tf, "cached": cached, "sig": sig, "cmd": None if cached else ["wsl", "-e", "bash", wsl(os.path.join(W, "run_%s.sh" % key))]}
print("trees:", {k: (len(v["aln"]), len(next(iter(v["aln"].values()))), "cached" if v["cached"] else "run") for k, v in jobs.items()}, flush=True)
t0 = time.time(); queue = [k for k, v in jobs.items() if v["cmd"]]; running = {}
queue.sort(key=lambda k: -len(jobs[k]["aln"]) * len(next(iter(jobs[k]["aln"].values()))))
while queue or running:
    while queue and len(running) < PAR:
        k = queue.pop(0); running[k] = subprocess.Popen(jobs[k]["cmd"])
    for k in list(running):
        if running[k].poll() is not None:
            if running[k].returncode != 0: raise RuntimeError("IQ-TREE failed for " + k)
            open(os.path.join(W, k + ".sha"), "w").write(jobs[k]["sig"]); print("  %s done after %.0f min" % (k, (time.time() - t0) / 60), flush=True); del running[k]
    time.sleep(5)
trees = {}
for key, j in jobs.items():
    nw = open(j["tf"], encoding="utf-8").read().strip(); leaves = re.findall(r"[(,]\s*([A-Z]\d+)\s*:", nw)
    assert sorted(leaves) == sorted(j["aln"]), "%s: tree leaves differ from the alignment rows" % key
    info = {"newick": nw, "n_taxa": len(j["aln"]), "n_columns": len(next(iter(j["aln"].values()))), "cached": j["cached"], "model": "LG+G4", "label": "%s · %s" % (LABEL[key.split("_")[1]], "structure-based" if key.startswith("struct") and meta["frames"][0]["kind"] == "structure" else "sequence-based")}
    p = os.path.join(W, key + ".iqtree")
    if os.path.exists(p):
        t = open(p, encoding="utf-8", errors="replace").read(); m = re.search(r"Log-likelihood of the tree: (-?[\d.]+)", t); info["logL"] = float(m.group(1)) if m else None
    trees[key] = info; print("%-12s %3d taxa %5d cols logL %s" % (key, info["n_taxa"], info["n_columns"], info.get("logL")), flush=True)
out = {"method": "IQ-TREE 2.3.6, fixed LG+G4, --fast search, SH-aLRT support (1000 replicates; values are SH-aLRT, not ultrafast bootstrap), seed %d; full-toxin trees use the columns occupied by >= 20 %% of the entries; domain trees use the columns of the module, entries with >= 100 (DI, DII) or 70 (DIII) residues of it" % SEED,
       "rooting": "stored unrooted; the page roots at the midpoint or on a chosen entry", "frames": {f["id"]: f["label"] for f in meta["frames"]}, "sets": {"A": "Cry toxin", **{k: LABEL[k] for k in SETS if k != "A"}}, "trees": trees}
json.dump(out, open(os.path.join(D, "trees.json"), "w", encoding="utf-8"), separators=(",", ":")); print("written trees.json", round(os.path.getsize(os.path.join(D, "trees.json")) / 1024), "KB")
