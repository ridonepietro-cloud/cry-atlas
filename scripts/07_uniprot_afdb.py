"""BPPRC sequence -> UniParc (MD5 of the identical sequence) -> UniProtKB accessions -> AlphaFold DB prediction (accepted only when the AFDB sequence equals the BPPRC sequence).
Output work/afdb/map.json {md5: {names, uniprot:[acc], reviewed:{acc:bool}, afdb:{acc, entry, version, pdbUrl, cifUrl, plddt} | None}} ; resumable.   Only metadata is fetched here (no model files)."""
import os, sys, re, json, hashlib, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
rows_seed = jload(os.path.join(WORK, "rows_seed.json"))                      # one row per distinct sequence (Cry3Aa10 resolved to the NCBI-identical record)
by = collections.OrderedDict()
for r in rows_seed:
    by.setdefault(hashlib.md5(r["seq"].encode()).hexdigest().upper(), {"seq": r["seq"], "names": []})["names"].extend(r["names"])
os.makedirs(os.path.join(WORK, "afdb"), exist_ok=True); OUTP = os.path.join(WORK, "afdb", "map.json"); out = jload(OUTP, {})
print("distinct sequences", len(by), "done", len(out), flush=True)
k = 0
for md5, d in by.items():
    if md5 in out: continue
    s = d["seq"]; accs = []; rev = {}
    if "X" in s:                                         # UniParc checksums are of the exact string: still try
        pass
    t = get("https://rest.uniprot.org/uniparc/search?query=%s&fields=upi,accession&size=5&format=json" % urllib.parse.quote("checksum:" + md5), key="uniprot", gap=0.25)
    for x in (json.loads(t or "{}").get("results") or []):
        for a in x.get("uniProtKBAccessions", []) or []:
            a = re.sub(r"\.\d+$", "", a)
            if a not in accs: accs.append(a)
    af = None
    for a in accs[:8]:
        j = get("https://alphafold.ebi.ac.uk/api/prediction/" + a, key="afdb", gap=0.15)
        p = (json.loads(j) or [None])[0] if j else None
        if p and p.get("pdbUrl") and p.get("uniprotSequence") == s:
            af = {"acc": a, "entry": p.get("entryId"), "version": p.get("latestVersion"), "pdbUrl": p.get("pdbUrl"), "cifUrl": p.get("cifUrl"), "plddt": p.get("globalMetricValue")}; break
    out[md5] = {"names": d["names"], "len": len(s), "uniprot": accs, "afdb": af}
    k += 1
    if k % 40 == 0:
        jdump(out, OUTP, indent=0); print("  done %d of %d; with UniProt %d, with AFDB model %d" % (len(out), len(by), sum(1 for v in out.values() if v["uniprot"]), sum(1 for v in out.values() if v["afdb"])), flush=True)
jdump(out, OUTP, indent=0)
print("finished: distinct", len(out), "with UniProt", sum(1 for v in out.values() if v["uniprot"]), "with AFDB model", sum(1 for v in out.values() if v["afdb"]))
