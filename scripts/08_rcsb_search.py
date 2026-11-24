"""RCSB sequence search (polymer entities >= 80 % identity, E <= 1e-20) for every distinct BPPRC Cry sequence -> work/pdb/rcsb_hits.json {md5: [[entity_id, identity, query_cov, subject_len, align_len, evalue]]}; resumable"""
import os, sys, re, json, hashlib, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
seqs = jload(os.path.join(DATA, "cry_seq_clean.json"))["seqs"]
by = collections.OrderedDict()
for n, s in seqs.items():
    by.setdefault(hashlib.md5(s.encode()).hexdigest().upper(), {"seq": s, "names": []})["names"].append(n)
os.makedirs(os.path.join(WORK, "pdb"), exist_ok=True); OUTP = os.path.join(WORK, "pdb", "rcsb_hits.json"); out = jload(OUTP, {})
print("distinct", len(by), "done", len(out), flush=True)
err = []; k = 0
for md5, d in by.items():
    if md5 in out: continue
    s = d["seq"].replace("X", "A")
    q = {"query": {"type": "terminal", "service": "sequence", "parameters": {"evalue_cutoff": 1e-20, "identity_cutoff": 0.8, "sequence_type": "protein", "value": s}},
         "return_type": "polymer_entity", "request_options": {"scoring_strategy": "sequence", "paginate": {"start": 0, "rows": 200}, "results_verbosity": "verbose"}}
    try:
        t = get("https://search.rcsb.org/rcsbsearch/v2/query", data=json.dumps(q).encode(), headers={"Content-Type": "application/json"}, key="rcsb", gap=0.6, tries=3)
    except Exception as e:
        print("  error for", d["names"][0], str(e)[:80], flush=True); err.append(md5); t = None
    hits = []
    if t:
        for h in json.loads(t).get("result_set", []):
            mc = (h.get("services") or [{}])[0].get("nodes", [{}])[0].get("match_context", [{}])[0]
            hits.append([h["identifier"], mc.get("sequence_identity"), mc.get("query_length"), mc.get("subject_length"), mc.get("alignment_length"), mc.get("evalue")])
    if md5 in err: continue
    out[md5] = hits; k += 1
    if k % 10 == 0:
        jdump(out, OUTP, indent=0); print("  done %d of %d; with hits %d" % (len(out), len(by), sum(1 for v in out.values() if v)), flush=True)
jdump(out, OUTP, indent=0)
print("finished; with hits", sum(1 for v in out.values() if v), "of", len(out))
