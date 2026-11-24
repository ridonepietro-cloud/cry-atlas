"""Targets of the genome-map step (work/cry/gc_targets.json): per row the protein accessions of its BPPRC entries (INSDC / RefSeq / UniProt-style), the UniProt accessions of identical sequences and the species.
Second pass (--locus): rows without an identical-protein record get an explicit location from the GenPept 'coded_by' of their entries (nucleotide accession + CDS range) -> work/cry/gc_targets_locus.json.
Run: py -I 18a_gc_targets.py [--locus]"""
import os, sys, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
rows = jload(os.path.join(WORK, "rows.json")); gp = jload(os.path.join(WORK, "ncbi", "genpept.json"))
os.makedirs(os.path.join(WORK, "cry"), exist_ok=True)
PROT = re.compile(r"^[A-Z]{3}\d{5}(\.\d+)?$|^[A-Z]{2}_\d+(\.\d+)?$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$")
if "--locus" not in sys.argv:
    t = {}
    for r in rows:
        accs = []
        for n, a in r["accessions"]:
            if PROT.match(a) and a not in accs: accs.append(a)
        for u in r.get("uniprot", [])[:3]:
            if u not in accs: accs.append(u)
        t[r["id"]] = {"accs": accs, "org": r["species"] or "Bacillus thuringiensis", "label": r["name"]}
    jdump(t, os.path.join(WORK, "cry", "gc_targets.json"), indent=0)
    print("targets", len(t), "with accessions", sum(1 for v in t.values() if v["accs"]))
else:
    done = jload(os.path.join(SITE, "data", "cry", "gc.json"), {"entries": {}})["entries"]
    out = {}
    for r in rows:
        if (done.get(r["id"]) or {}).get("loci"): continue
        for n, cb in r.get("coded_by", []):
            m = re.match(r"^(?:complement\()?(?:join\()?([A-Z]{1,2}_?\d+(?:\.\d+)?):<?(\d+)\.\.>?(\d+)", cb)
            if m:
                strand = "-" if cb.startswith("complement") else "+"
                out[r["id"]] = {"accs": [], "org": r["species"] or "Bacillus thuringiensis", "label": r["name"], "locus": [m.group(1), int(m.group(2)), int(m.group(3)), strand, "CDS location in the GenPept record of %s (coded_by)" % n]}; break
    jdump(out, os.path.join(WORK, "cry", "gc_targets_locus.json"), indent=0); print("rows with a locus from coded_by:", len(out))
