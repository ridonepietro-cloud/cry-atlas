"""Parse the BPPRC structures page (raw/web2/structures_.html): the table has malformed rows (only the first <tr> opens), so cells are read in groups of 8.
Output: data/bpprc_cry_structures.tsv (name(s), old name, accession(s), pdb, modified, comment, year, doi) for Cry-named entries"""
import os, sys, re, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
t = open(os.path.join(CRY, "raw", "web2", "structures_.html"), encoding="utf8", errors="replace").read()
i = t.find("<tbody"); body = t[i:t.find("</tbody>", i)]
tds = re.findall(r"<TD[^>]*>(.*?)</TD>", body, re.S | re.I)
cl = [re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", x)).strip() for x in tds]
print(len(cl), "cells ->", len(cl) / 8.0, "rows")
rows = [cl[k:k + 8] for k in range(0, len(cl) - 7, 8)]
allr = rows
cry = [r for r in rows if re.search(r"\bCry\d|Cry [Cc]him|Cry [Vv]ar|^Cry\b|cry", r[0], re.I)]
print("all rows", len(rows), "Cry-like rows", len(cry))
for r in cry: print(r)
with open(os.path.join(DATA, "bpprc_cry_structures.tsv"), "w", encoding="utf8", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["name", "old_name", "accession", "pdb", "modified", "comment", "year", "doi"]); w.writerows(cry)
print("non-Cry classes in file:", collections.Counter(re.match(r"[A-Za-z]+", r[0]).group(0)[:3] for r in rows if r not in cry))
