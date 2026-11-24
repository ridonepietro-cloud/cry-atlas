"""Parse the BPPRC Cry category page (raw/web/cry_cat.html) into cry_atlas/data/bpprc_cry_table.tsv
Rows with a checkbox = in the downloadable set; rows with '*' in the first cell = no checkbox (named, not in the cart/download set)."""
import re, csv, os, collections
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
t = open(os.path.join(ROOT, "raw", "web", "cry_cat.html"), encoding="utf8", errors="replace").read()
title = re.search(r"<h5>(.*?)</h5>", t, re.S).group(1).strip()
body = t[t.find("<tbody>"):t.find("</tbody>")]
rows = []
for tr in re.findall(r"<tr>(.*?)</tr>", body, re.S):
    tds = re.findall(r"<TD[^>]*>(.*?)(?:</TD>|</td>)", tr, re.S | re.I)
    cl = [re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", x)).strip() for x in tds]
    cart = 'name="name"' in tr
    nm = re.search(r'value="([^"]+)"', tr).group(1) if cart else cl[1]
    accm = re.search(r'/protein/([^"]+)"', tr)
    rows.append({"name": nm.strip(), "in_cart": 1 if cart else 0, "old_name": cl[2], "other_name": cl[3], "accession": accm.group(1).strip() if accm else cl[4], "year": cl[5]})
print(title, len(rows), "in cart:", sum(r["in_cart"] for r in rows))
os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
with open(os.path.join(ROOT, "data", "bpprc_cry_table.tsv"), "w", encoding="utf8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
c = collections.Counter(re.match(r"Cry(\d+)", r["name"]).group(1) if re.match(r"Cry(\d+)", r["name"]) else "other" for r in rows)
print("families:", len(c), sorted(c.items(), key=lambda x: (not x[0].isdigit(), int(x[0]) if x[0].isdigit() else 0))[:80])
acc = collections.Counter()
for r in rows:
    a = r["accession"]
    k = "ncbi-like" if re.fullmatch(r"[A-Z]{1,3}_?\d+(\.\d+)?", a) else ("patent/other" if a else "none")
    acc[k] += 1
print(acc)
print([r for r in rows if not re.fullmatch(r"[A-Z]{1,3}_?\d+(\.\d+)?", r["accession"])][:15])
