"""NCBI GenPept of every BPPRC accession (organism, strain, plasmid, protein name, coded_by nucleotide record, NCBI sequence) -> work/ncbi/genpept.json ; raw text in RAW/ncbi_gp/
Accessions that are not protein records (WGS contigs, patents, UniProt-only) are listed as 'no record'.   Run: py -I 04_ncbi_genpept.py"""
import os, sys, re, json, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
tab = list(csv.DictReader(open(os.path.join(DATA, "bpprc_cry_table.tsv"), encoding="utf8"), delimiter="\t"))
os.makedirs(os.path.join(WORK, "ncbi"), exist_ok=True); os.makedirs(os.path.join(RAW, "ncbi_gp"), exist_ok=True)
OUTP = os.path.join(WORK, "ncbi", "genpept.json"); out = jload(OUTP, {})
accs = sorted({r["accession"] for r in tab if re.fullmatch(r"[A-Z]{1,3}\d{5,}(\.\d+)?|[A-Z]{2}_\d+(\.\d+)?|[OPQ]\d[A-Z0-9]{3}\d", r["accession"])})
todo = [a for a in accs if a not in out]
print("accessions", len(accs), "to fetch", len(todo), flush=True)


def parse_gp(txt):
    recs = {}
    for block in txt.split("\n//")[:-1] if "\n//" in txt else []:
        m = re.search(r"^VERSION\s+(\S+)", block, re.M)
        if not m: continue
        ver = m.group(1); acc = re.search(r"^ACCESSION\s+(\S+)", block, re.M).group(1)
        r = {"version": ver, "accession": acc}
        r["definition"] = re.sub(r"\s+", " ", re.search(r"^DEFINITION\s+(.*?)\n(?=\S)", block, re.S | re.M).group(1)).strip()
        so = re.search(r"^SOURCE\s+(.*?)\n\s+ORGANISM\s+(.*?)\n(?=\S)", block, re.S | re.M)
        if so:
            r["source"] = re.sub(r"\s+", " ", so.group(1)).strip(); lin = re.sub(r"\s+", " ", so.group(2)); r["organism"] = lin.split(" Bacteria;")[0].strip() if " Bacteria;" in lin else lin.split(";")[0].strip(); r["lineage"] = lin
        lm = re.search(r"^LOCUS\s+\S+\s+(\d+)\s*aa", block, re.M)
        r["length"] = int(lm.group(1)) if lm else None
        for key in ("strain", "plasmid", "isolate", "serovar", "sub_species", "country", "note", "host", "db_xref", "coded_by"):
            vals = re.findall(r"/%s=\"([^\"]*)\"" % key, block, re.S)
            if vals: r[key] = [re.sub(r"\s+", " ", v) for v in vals]
        cb = re.search(r"/coded_by=\"([^\"]*)\"", block, re.S)
        if cb: r["coded_by"] = re.sub(r"\s+", "", cb.group(1))
        og = re.search(r"^ORIGIN\s*\n(.*)", block, re.S | re.M)
        r["seq"] = re.sub(r"[^a-z]", "", og.group(1)).upper() if og else ""
        r["dbsource"] = (re.search(r"^DBSOURCE\s+(.*)", block, re.M) or [None, ""])[1] if re.search(r"^DBSOURCE", block, re.M) else ""
        recs[ver] = r; recs[acc] = r
    return recs


for i in range(0, len(todo), 100):
    chunk = todo[i:i + 100]
    body = urllib.parse.urlencode({"db": "protein", "id": ",".join(chunk), "rettype": "gp", "retmode": "text", "tool": "atlas-cry-pipeline"}).encode()
    t = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", data=body, gap=0.4)
    open(os.path.join(RAW, "ncbi_gp", "gp_%04d.txt" % i), "w", encoding="utf-8").write(t or "")
    recs = parse_gp(t or "")
    for a in chunk:
        base = a.split(".")[0]
        r = recs.get(a) or recs.get(base)
        out[a] = r if r else {"missing": True}
    print("batch", i, "records parsed", len(set(id(v) for v in recs.values())), "of", len(chunk), flush=True)
    jdump(out, OUTP, indent=0)
jdump(out, OUTP, indent=0)
print("found", sum(1 for v in out.values() if not v.get("missing")), "missing", [a for a, v in out.items() if v.get("missing")])
