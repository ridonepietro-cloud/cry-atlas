"""Entity sequences of the PDB candidates (RCSB GraphQL, 100 entries per request) and local matching to the Cry rows.
A candidate entity is KEPT when it is >= 90 % identical to a row over >= 70 % of the entity length (5-mer prefilter, then Biopython global alignment with free end gaps).
Output work/pdb/entities.json (all candidate entities) and work/pdb/matches.json {pdb_id: [{entity, chains, desc, organism, len, rows: [[row_id, identity, covered_fraction]]}]}.   Run: py -I 08c_pdb_entities.py"""
import os, sys, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
cand = jload(os.path.join(WORK, "pdb", "candidates.json")); ids = sorted(cand)
print("candidate entries", len(ids), flush=True)
EP = os.path.join(WORK, "pdb", "entities.json"); ent = jload(EP, {})
Q = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title} rcsb_accession_info{initial_release_date} exptl{method} polymer_entities{rcsb_id rcsb_polymer_entity{pdbx_description} entity_poly{pdbx_seq_one_letter_code_can} rcsb_polymer_entity_container_identifiers{auth_asym_ids reference_sequence_identifiers{database_name database_accession}} rcsb_entity_source_organism{ncbi_scientific_name}}}}"""
todo = [i for i in ids if i not in ent]
for k in range(0, len(todo), 100):
    chunk = todo[k:k + 100]
    t = get("https://data.rcsb.org/graphql", data=json.dumps({"query": Q, "variables": {"ids": chunk}}).encode(), headers={"Content-Type": "application/json"}, key="rcsbdata", gap=0.4, tries=3)
    d = json.loads(t)["data"]["entries"] if t else []
    for e in d or []:
        ent[e["rcsb_id"]] = {"title": e["struct"]["title"], "date": e["rcsb_accession_info"]["initial_release_date"], "method": [x["method"] for x in e["exptl"]],
                             "entities": [{"id": p["rcsb_id"], "desc": (p["rcsb_polymer_entity"] or {}).get("pdbx_description"), "seq": (p["entity_poly"] or {}).get("pdbx_seq_one_letter_code_can", "").replace("\n", ""),
                                           "chains": (p["rcsb_polymer_entity_container_identifiers"] or {}).get("auth_asym_ids") or [],
                                           "uniprot": [r["database_accession"] for r in ((p["rcsb_polymer_entity_container_identifiers"] or {}).get("reference_sequence_identifiers") or []) if r["database_name"] == "UniProt"],
                                           "organism": [o["ncbi_scientific_name"] for o in (p["rcsb_entity_source_organism"] or [])]} for p in (e["polymer_entities"] or [])]}
    print("  entities for", min(k + 100, len(todo)), "of", len(todo), flush=True)
    jdump(ent, EP)
jdump(ent, EP)
rows = jload(os.path.join(WORK, "rows_seed.json"))
from Bio.Align import PairwiseAligner
al = PairwiseAligner(); al.mode = "global"; al.match_score = 2; al.mismatch_score = -1; al.open_gap_score = -6; al.extend_gap_score = -1; al.end_gap_score = 0


def kset(s, k=5): return {s[i:i + k] for i in range(len(s) - k + 1)}


rk = {r["id"]: kset(r["seq"]) for r in rows}; rs = {r["id"]: r["seq"] for r in rows}
out = {}
for pid, e in ent.items():
    for en in e["entities"]:
        s = en["seq"]
        if len(s) < 100: continue
        k = kset(s); best = []
        for rid, ks in rk.items():
            sh = len(k & ks) / max(1, min(len(k), len(ks)))
            if sh >= 0.25: best.append((sh, rid))
        best.sort(reverse=True); res = []
        for sh, rid in best[:60]:
            a = al.align(rs[rid], s)[0]; same = cov = 0
            for (r0, r1), (c0, c1) in zip(*a.aligned):
                cov += r1 - r0; same += sum(1 for q in range(r1 - r0) if rs[rid][r0 + q] == s[c0 + q])
            ident = same / cov if cov else 0; frac = cov / len(s)
            if ident >= 0.90 and frac >= 0.70: res.append([rid, round(ident, 4), round(frac, 3)])
        if res:
            out.setdefault(pid, []).append({"entity": en["id"], "chains": en["chains"], "desc": en["desc"], "organism": en["organism"], "len": len(s), "uniprot": en["uniprot"], "title": e["title"], "date": e["date"], "method": e["method"], "rows": sorted(res, key=lambda x: -x[1])})
jdump(out, os.path.join(WORK, "pdb", "matches.json"), indent=0)
print("entries with a Cry-row match:", len(out), "of", len(ent))
for pid, v in sorted(out.items()):
    for m in v: print(" ", pid, m["entity"], m["desc"], "len", m["len"], "rows", len(m["rows"]), "best", m["rows"][0], "| " + ";".join(m["organism"])[:40])
