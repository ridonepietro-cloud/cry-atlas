"""PDB candidates for the Cry rows without the (server-side timing out) sequence search:
 (1) RCSB attribute search: polymer entities whose UniProt accession is one of the UniProt accessions of the rows (work/afdb/map.json), in chunks of 200 accessions;
 (2) RCSB full-text search for Cry-toxin terms (entries), filtered later by sequence identity to the rows;
 (3) the entries of the BPPRC structure page.
Output work/pdb/candidates.json {pdb_id: [sources]}.    Run: py -I 08b_rcsb_candidates.py"""
import os, sys, re, json, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
os.makedirs(os.path.join(WORK, "pdb"), exist_ok=True)
mp = jload(os.path.join(WORK, "afdb", "map.json"))
accs = sorted({a for v in mp.values() for a in v.get("uniprot", [])})
print("UniProt accessions", len(accs), flush=True)
cand = {}


def add(i, src): cand.setdefault(i.upper(), set()).add(src)


def search(q, rows=1000):
    body = {"query": q, "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": rows}, "results_verbosity": "compact"}}
    t = get("https://search.rcsb.org/rcsbsearch/v2/query", data=json.dumps(body).encode(), headers={"Content-Type": "application/json"}, key="rcsb", gap=0.5, tries=3)
    return json.loads(t).get("result_set", []) if t else []


for k in range(0, len(accs), 150):
    chunk = accs[k:k + 150]
    q = {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "operator": "in", "value": chunk}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name", "operator": "exact_match", "value": "UniProt"}}]}
    try:
        for r in search(q): add(r if isinstance(r, str) else r["identifier"], "uniprot")
    except Exception as e:
        print("  uniprot chunk", k, "error", str(e)[:100], flush=True)
print("by UniProt accession:", len(cand), flush=True)
for term in ("delta-endotoxin", "pesticidal crystal protein", "insecticidal crystal protein", "Cry toxin Bacillus thuringiensis", "Cry1Ac", "Cry3Aa", "Cry2Aa", "Cry4Ba", "Cry5B", "Cry8", "Cry11Aa", "Cry1Ca", "Cry1Fa", "Cry1Ba", "Cry7", "Cry34", "protoxin crystal Bacillus thuringiensis"):
    try:
        for r in search({"type": "terminal", "service": "full_text", "parameters": {"value": term}}, rows=300): add(r if isinstance(r, str) else r["identifier"], "text:" + term)
    except Exception as e:
        print("  text", term, "error", str(e)[:100], flush=True)
print("with full-text hits:", len(cand), flush=True)
for r in csv.DictReader(open(os.path.join(DATA, "bpprc_cry_structures.tsv"), encoding="utf-8"), delimiter="\t"):
    for x in re.findall(r"\b[0-9][A-Za-z0-9]{3}\b", r["pdb"]): add(x, "bpprc")
print("with BPPRC list:", len(cand))
jdump({k: sorted(v) for k, v in sorted(cand.items())}, os.path.join(WORK, "pdb", "candidates.json"), indent=0)
