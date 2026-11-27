"""Ligand-site conservation (2026-11-27): conservation of the ligand-contacting residues of a Cry structure with a bound sugar.
Demonstration with PDB 4ARY (Cry1Ac with N-acetylgalactosamine, from the BPPRC structures table): sugar atoms -> protein residues within CUTOFF -> row index of the matching Atlas row -> column of the structure-based frame ->
conservation of that column over all rows, over the nematode-active rows and over the Lepidoptera-active rows.   Reusable for another structure: py -I 23_ligand_site.py <mmCIF path> <row id> [ligand ids comma separated]
Output: printed table + docs/ligand_site_<pdb>.tsv.     Needs RAW/pdb/4ARY.cif (12_fetch_structures.py)."""
import os, sys, json, math, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
CIF = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RAW, "pdb", "4ARY.cif"); ROWID = sys.argv[2] if len(sys.argv) > 2 else None
LIG = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else None; CUTOFF = 4.5
rows = jload(os.path.join(SITE, "data", "cry", "rows.json")); meta = jload(os.path.join(SITE, "data", "cry", "meta.json")); R = {r["id"]: r for r in rows}
src = {r["id"]: r for r in jload(os.path.join(WORK, "rows.json"))}
pdb = os.path.splitext(os.path.basename(CIF))[0].upper()
if ROWID is None:
    mt = jload(os.path.join(WORK, "pdb", "matches.json"))[pdb]; ROWID = mt[0]["rows"][0][0]
print("structure", pdb, "-> Atlas row", ROWID, R[ROWID]["label"], "| sequence length", len(src[ROWID]["seq"]))
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
d = MMCIF2Dict(CIF); g = lambda k: d.get(k, []) if isinstance(d.get(k, []), list) else [d.get(k)]
grp, comp, chain, seqid, ax, ay, az, aname = (g("_atom_site." + k) for k in ("group_PDB", "label_comp_id", "auth_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "label_atom_id"))
het = collections.Counter((c, ch) for gp, c, ch in zip(grp, comp, chain) if gp == "HETATM" and c != "HOH")
print("non-water HETATM groups (comp, chain): count", dict(het))
sugar = LIG or {c for (c, ch) in het if c in ("NGA", "A2G", "GLA", "GAL", "NAG", "GM0", "GLC", "MAN", "BGC", "GAD", "A2G")}
print("ligand ids used:", sorted(sugar))
lig_atoms = [(float(ax[i]), float(ay[i]), float(az[i]), comp[i], chain[i], seqid[i]) for i in range(len(grp)) if grp[i] == "HETATM" and comp[i] in sugar]
print("ligand atoms:", len(lig_atoms), "| ligand residues:", sorted({(a[3], a[4], a[5]) for a in lig_atoms}))
contacts = collections.defaultdict(lambda: [99.0, set()])
for i in range(len(grp)):
    if grp[i] != "ATOM": continue
    p = (float(ax[i]), float(ay[i]), float(az[i]))
    for (x, y, z, lc, lch, ls) in lig_atoms:
        dd = math.sqrt((p[0] - x) ** 2 + (p[1] - y) ** 2 + (p[2] - z) ** 2)
        if dd <= CUTOFF:
            k = (chain[i], int(seqid[i]), comp[i]); c = contacts[k]; c[0] = min(c[0], dd); c[1].add(aname[i])
print("protein residues within %.1f A of the ligand: %d" % (CUTOFF, len(contacts)))
# map to the row sequence
sch = CC.cif_scheme(CIF); seq = src[ROWID]["seq"]; chains_used = sorted({k[0] for k in contacts}); mapping = {}
for ch in chains_used:
    res = CC.read_chains(CIF, {ch})[ch]; idx, ident, mode = CC.place_chain(res, seq, ch, sch)
    print("chain %s: %d residues placed on the row sequence, identity %.3f (%s)" % (ch, sum(1 for i in idx if i >= 0), ident, mode))
    for x, i in zip(res, idx):
        if i >= 0: mapping[(ch, x["resnum"])] = i
frame = jload(os.path.join(SITE, "data", "cry", meta["frames"][0]["file"]))["seqs"]; cols = CC.col_index(frame[ROWID])
import base64
modname = meta["module_ids"]; modcode = base64.b64decode(jload(os.path.join(SITE, "data", "cry", "residues.json"))[ROWID]["A"]["mod"])
spec = jload(os.path.join(WORK, "specificity_breadth.json")); nem = set(spec["nematode_rows"])
sp = jload(os.path.join(WORK, "specificity.json")); lep = {r["id"] for r in rows if any(n in sp and "Lepidoptera" in sp[n]["yes"] for n in (r.get("extra") or {}).get("names", []))}
cry1 = {r["id"] for r in rows if (r.get("extra") or {}).get("family") == "Cry1"}
print("nematode-active rows %d | Lepidoptera-active rows %d | Cry1 rows %d | all rows %d" % (len(nem), len(lep), len(cry1), len(rows)))
ref_aa = {}
def column_stats(c, refaa, subset=None):
    chars = [frame[r["id"]][c - 1] for r in rows if (subset is None or r["id"] in subset)]
    res = [x for x in chars if x != "-"]
    if not res: return None
    same = sum(1 for x in res if x == refaa) / float(len(res)); cnt = collections.Counter(res)
    ent = -sum((n / float(len(res))) * math.log(n / float(len(res)), 2) for n in cnt.values())
    return {"n": len(res), "of": len(chars), "same_as_ref": same, "top": cnt.most_common(3), "entropy": ent}
out = []
print("%-8s %-4s %-6s %-5s | %-34s | all rows: ref-aa conserved (n)  top residues                      | nematode-active: conserved (n) | Lepidoptera-active: conserved (n)" % ("residue", "aa", "module", "d(A)", "atoms"))
merged = {}                                                                 # one entry per residue number: the sugar is bound in several chains of the crystal; shortest distance, atoms of all copies
for (ch, rn, aa3), (dmin, atoms) in contacts.items():
    m = merged.setdefault((rn, aa3), [99.0, set(), None]); m[0] = min(m[0], dmin); m[1] |= atoms
    if m[2] is None and (ch, rn) in mapping: m[2] = mapping[(ch, rn)]
print("distinct contact residues (merged over the chains): %d" % len(merged))
for (rn, aa3), (dmin, atoms, i) in sorted(merged.items()):
    ch = '*'
    if i is None: print("residue", rn, aa3, "not placed"); continue
    c = cols[i]; one = seq[i]; mod = modname.get(str(modcode[i]), "?")
    a = column_stats(c, one); b = column_stats(c, one, nem); l = column_stats(c, one, lep); k1 = column_stats(c, one, cry1)
    f = lambda s: "%.2f (%d)" % (s["same_as_ref"], s["n"]) if s else "no residues"
    print("%-8s %-4s %-6s %-5.1f | %-34s | %-22s %-34s | %-30s | %s" % ("%s%d" % (one, rn), aa3, mod, dmin, ",".join(sorted(atoms))[:34], f(a), ",".join("%s%d" % t for t in a["top"]) if a else "", f(b), f(l) + "  | Cry1: " + f(k1)))
    out.append([one + str(rn), aa3, mod, round(dmin, 2), c, i + 1, a["same_as_ref"] if a else "", a["n"] if a else "", b["same_as_ref"] if b else "", b["n"] if b else "", l["same_as_ref"] if l else "", l["n"] if l else "", ";".join("%s%d" % t for t in a["top"]) if a else ""])
# context: conservation of all columns of the same module family (median over residues of the reference row)
allc = [column_stats(cols[j], seq[j]) for j in range(len(seq))]
vals = sorted(s["same_as_ref"] for s in allc if s); print("context: all %d residues of %s: median fraction of rows with the same residue as the reference = %.2f; the ligand-contact residues: median %.2f" % (len(vals), ROWID, vals[len(vals) // 2], sorted(o[6] for o in out if o[6] != "")[len([o for o in out if o[6] != ""]) // 2] if out else float("nan")))
import csv
p = os.path.join(CRY, "docs", "ligand_site_%s.tsv" % pdb.lower())
with open(p, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["residue", "comp", "module", "min_dist_A", "frame_col", "row_index_1based", "frac_same_all", "n_all", "frac_same_nematode_active", "n_nematode", "frac_same_lepidoptera_active", "n_lepidoptera", "top_residues_all_rows"]); w.writerows(out)
print("written", p)
