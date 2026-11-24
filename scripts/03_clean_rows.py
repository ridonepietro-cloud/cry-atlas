"""Clean BPPRC sequences, merge exact duplicates -> data/cry_seq_clean.json (per BPPRC name: cleaned seq, issues) and data/cry_rows_seed.json (one row per distinct sequence)."""
import os, re, csv, json, collections, hashlib
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = os.path.join(ROOT, "raw", "fasta", "cry_full.fasta")
recs = []; k = None
for l in open(F, encoding="utf8"):
    l = l.rstrip("\n")
    if l.startswith(">"): k = [l[1:].strip(), []]; recs.append(k)
    elif k is not None: k[1].append(l)
names = [r[0] for r in recs]
dup = [n for n, c in collections.Counter(names).items() if c > 1]
print("records", len(recs), "duplicate names", dup)
for n in dup:
    ss = ["".join(r[1]) for r in recs if r[0] == n]; print(n, [len(s) for s in ss], "identical" if len(set(ss)) == 1 else "DIFFERENT")
tab = {r["name"]: r for r in csv.DictReader(open(os.path.join(ROOT, "data", "bpprc_cry_table.tsv"), encoding="utf8"), delimiter="\t")}
hol = set(l[1:].strip() for l in open(os.path.join(ROOT, "raw", "fasta", "holotype.fasta"), encoding="utf8") if l.startswith(">"))
clean = collections.OrderedDict(); issues = {}
for n, parts in recs:
    raw = "\n".join(parts)
    s = re.sub(r"\s+", "", raw)
    iss = []
    if re.search(r"[a-z]", s): iss.append("lowercase letters in the downloaded sequence (uppercased)")
    if re.search(r"[ \t]", "".join(parts)): iss.append("blank inside a sequence line (removed)")
    s = s.upper()
    if re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", s): iss.append("X / non-standard residue(s): " + "".join(sorted(set(re.findall(r"[^ACDEFGHIKLMNPQRSTVWY]", s)))))
    if n in clean:
        if clean[n] != s: iss.append("second record of the same name with a different sequence")
        continue
    clean[n] = s
    if iss: issues[n] = iss
print("issues:", len(issues)); 
for n, v in list(issues.items())[:40]: print(" ", n, v)
# the lowercase one
print(recs[[r[0] for r in recs].index("Cry8Ba2")][1][:3])
# exact duplicates
byseq = collections.defaultdict(list)
for n, s in clean.items(): byseq[s].append(n)
print("distinct sequences", len(byseq), "groups with >1 name", sum(1 for v in byseq.values() if len(v) > 1), "names in them", sum(len(v) for v in byseq.values() if len(v) > 1))
print("largest groups", sorted(((len(v), v[:4]) for v in byseq.values()), reverse=True)[:8])
json.dump({"seqs": clean, "issues": issues}, open(os.path.join(ROOT, "data", "cry_seq_clean.json"), "w"), indent=0)
