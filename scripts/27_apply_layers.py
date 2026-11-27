"""Post-build layers of the Cry pack (run AFTER 19_build_cry.py, which rewrites rows / features / meta; idempotent):
  * row fields 'activity' (nematode-active / insect-active / ...) and 'domtree' (domain trees concordant / one domain differs / ...) + one info line each      -> rows.json
  * two colour-by / filter attributes 'Activity group' and 'Domain trees'                                                                                      -> features.json
  * marks (diamonds / stars on the alignment columns of every row that has a residue there): GalNAc-contacting residues of Cry1Ac (4ARY) and the residue columns shared by the nematode toxins -> meta.json
  * tanglegram default = Domain II against Domain III (domain swaps)                                                                                           -> meta.json
  * docs/nematode_signature_contacts.tsv, docs/nematode_vs_insect_separation.md, docs/domain_tree_discordance.md
Inputs: work/separation.json (24_groups_and_separation.py), work/domain_discordance.json (25_domain_trees.py), docs/ligand_site_4ary.tsv (23_ligand_site.py).     Run: py -I 27_apply_layers.py"""
import os, sys, csv, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
D = os.path.join(SITE, "data", "cry")
rows = jload(os.path.join(D, "rows.json")); feat = jload(os.path.join(D, "features.json")); meta = jload(os.path.join(D, "meta.json"))
sep = jload(os.path.join(WORK, "separation.json")); dd = jload(os.path.join(WORK, "domain_discordance.json"))
lig = list(csv.DictReader(open(os.path.join(CRY, "docs", "ligand_site_4ary.tsv"), encoding="utf-8"), delimiter="\t"))
R = {r["id"]: r for r in rows}
ACT_COL = {"Nematode-active": "#c51b7d", "Insect-active": "#2f8a57", "Other organism only": "#8a5cc4", "Tested, none active": "#8a94a0", "No data": "#b9bec4"}
DOM_COL = {"concordant": "#2f8a57", "Domain I differs": "#d49a14", "Domain II differs": "#3b7cbd", "Domain III differs": "#d64045", "all three differ": "#8a5cc4", "not assessable": "#b9bec4"}
PREF = ("Activity group (BPPRC", "Domain trees (")


def sp_list(d, n=4):
    return "; ".join("%s (%s%s)" % (o, ", ".join(s[:n]), " ..." if len(s) > n else "") for o, s in sorted(d.items(), key=lambda kv: -len(kv[1])))


for r in rows:
    i = r["id"]; g = sep["groups"].get(i, "No data"); r["activity"] = g
    c = dd["category"].get(i)
    r["domtree"] = "not in all domain trees" if c is None else ("not assessable" if c == "mixed neighbourhood" else c)
    r["info"] = [x for x in r.get("info", []) if not x.startswith(PREF)]
    ev = sep["evidence"].get(i)
    if ev:
        t = "Activity group (BPPRC specificity): %s. Active: %s." % (g, sp_list(ev["active_orders"]) or "none recorded")
        if ev["inactive_orders"]: t += " Inactive records: %s." % sp_list(ev["inactive_orders"], 3)
    else: t = "Activity group (BPPRC specificity): no specificity record for any of its names."
    r["info"].append(t)
    if c is not None:
        gr = dd["groups"][i]; cl = dd["closest"][i]; nm = {"I": "DI", "II": "DII", "III": "DIII"}
        txt = " / ".join("%s %s" % (d, gr[d][0]) for d in ("DI", "DII", "DIII"))
        if r["domtree"] == "not assessable": r["info"].append("Domain trees (nearest relatives): the 10 nearest relatives in at least one domain tree are not dominated by one named group, so a domain swap cannot be judged here.")
        else: r["info"].append("Domain trees (nearest relatives): group of the 10 nearest relatives in the DI / DII / DIII tree = %s: %s. Closest entry per tree: %s." % (txt, r["domtree"], " / ".join(R[cl[d]]["label"].split(" ")[0] for d in ("DI", "DII", "DIII"))))
json.dump(rows, open(os.path.join(D, "rows.json"), "w", encoding="utf-8"), separators=(",", ":"))
# ---------- features
cnt_a = collections.Counter(r["activity"] for r in rows); cnt_d = collections.Counter(r["domtree"] for r in rows)
def attr(key, label, short, text, cnt, col, order):
    vals = [v for v in order if cnt.get(v)] + [v for v in cnt if v not in order]
    return {"key": key, "label": label, "short": short, "text": text, "values": [{"id": v, "label": v, "colour": col.get(v, "#b9bec4"), "n": cnt[v]} for v in vals]}
A1 = attr("activity", "Activity group (BPPRC specificity)", "Activity", "Nematode-active = an active record against nematodes; Insect-active = an active record against an insect order and none against nematodes; no data = no specificity record for any name of the sequence", cnt_a, ACT_COL, list(ACT_COL))
cnt_d2 = collections.Counter(r["domtree"] for r in rows); DOM_COL["not in all domain trees"] = "#b9bec4"
A2 = attr("domtree", "Domain trees: does one domain differ?", "Domain trees", "Group (BPPRC secondary rank, e.g. Cry1B) of the 10 nearest relatives in the Domain I, II and III trees; 'Domain X differs' = the other two trees agree; not assessable = no dominant group among the nearest relatives", cnt_d2, DOM_COL, list(DOM_COL))
feat["attrs"] = [A1, A2] + [a for a in feat["attrs"] if a["key"] not in ("activity", "domtree")]
json.dump(feat, open(os.path.join(D, "features.json"), "w", encoding="utf-8"), separators=(",", ":"))
# ---------- marks and tanglegram
sig = sep["top_columns"]; sigcols = sorted({o["col"] for o in sig if o["freq_pos"] >= 0.8 and o["freq_neg"] <= 0.1})
ligcols = sorted({int(x["frame_col"]) for x in lig})
marks = [m for m in meta.get("marks", []) if m["id"] not in ("galNAc4ary", "nemsig")]
marks.append({"id": "galNAc4ary", "label": "GalNAc-contacting residues of Cry1Ac (PDB 4ARY)", "sub": "A", "cols": ligcols, "color": "#c2185b", "glyph": "diamond",
              "note": "GalNAc site of Cry1Ac: the %d residues within 4.5 A of the bound N-acetylgalactosamine in 4ARY (Domain III: %s), carried to every row by the alignment column (structure-based frame). Not conserved across the Atlas: %s." % (len(ligcols), ", ".join(x["residue"] for x in lig), "see docs/ligand_site_4ary.tsv"), "source": "scripts/23_ligand_site.py"})
marks.append({"id": "nemsig", "label": "Residue columns shared by the nematode toxins", "sub": "A", "cols": sigcols, "color": "#00897b", "glyph": "star",
              "note": "Nematode signature (candidate): %d alignment columns where >= 80 %% of the nematode-active toxins of the Cry5, 12, 13, 14, 21 and 95 families carry the same residue and <= 10 %% of the insect-active toxins do. Holds when a whole nematode family is left out (mean AUC %.2f against %.2f for random families, p = %.3f), but nematode-active toxins are close relatives, so a lineage marker is not shown to be functional." % (len(sigcols), sep["columns"]["auc_mean"], sep["columns"]["null_mean"], sep["columns"]["p"]), "source": "scripts/24_groups_and_separation.py"})
meta["marks"] = marks; meta["tangle"] = {"default": ["struct_DII", "struct_DIII"]}
json.dump(meta, open(os.path.join(D, "meta.json"), "w", encoding="utf-8"), separators=(",", ":"))
# ---------- docs
with open(os.path.join(CRY, "docs", "nematode_signature_contacts.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["frame_col_1", "frame_col_2", "module_1", "module_2", "freq_nematode", "freq_insect", "Cry5Ba1_residue_1", "Cry5Ba1_residue_2"])
    for o in sep["top_contacts"]: w.writerow([o["cols"][0], o["cols"][1], o["modules"][0], o["modules"][1], round(o["freq_pos"], 3), round(o["freq_neg"], 3), o["cry5Ba1_pos"][0], o["cry5Ba1_pos"][1]])
cs, ct = sep["columns"], sep["contacts"]; cc = dd["counts"]; assess = sum(v for k, v in cc.items() if k != "mixed neighbourhood")
top = [x for x in sorted(dd["category"], key=lambda x: (-({"all three differ": 2, "Domain I differs": 1, "Domain II differs": 1, "Domain III differs": 1}.get(dd["category"][x], 0)), R[x]["label"])) if dd["category"][x] not in ("concordant", "mixed neighbourhood")]
open(os.path.join(CRY, "docs", "nematode_vs_insect_separation.md"), "w", encoding="utf-8").write("""# Nematode-active versus insect-active Cry toxins: which columns and contacts separate them (2026-11-27)

Groups from the BPPRC specificity export (per sequence, union over its names): %s.
Positives = %d nematode-active toxins of six families (%s); negatives = %d insect-active toxins of %d families.

**Alignment columns** (structure-based frame, columns with >= 30 %% occupancy): %d residue-column features pass all criteria (>= 80 %% of the positives, <= 10 %% of the negatives, in the majority of the rows of >= 5 of the 6 nematode families).
Leave-one-nematode-family-out AUC per held-out family: %s, mean %.2f; random pseudo-nematode families (%d draws): mean %.2f, 95th percentile %.2f, p = %.3f.

**Contacts** (stored contacts of the primary models, column pairs present in >= 3 %% of the rows): %d column pairs pass the same criteria. Leave-one-family-out mean AUC %.2f; random families mean %.2f, 95th percentile %.2f, p = %.3f.
Most of the passing contacts lie in the C-terminal (protoxin) module: a hint that the nematode toxins share a protoxin architecture that the insect-active toxins in the sample do not share at the same aligned positions.

Caveats: 17 positives in 6 families (Cry12, Cry13 and Cry95 are single entries); the nematode toxins are close relatives, so a signature is first a lineage marker; the family-label permutation test controls for this only in part; the activity labels are BPPRC records, which are incomplete (461 of 670 sequences have no specificity record).
Files: nematode_signature_columns.tsv, nematode_signature_contacts.tsv; scripts 24_groups_and_separation.py, 27_apply_layers.py.
""" % (dict(sep["group_counts"]), len(sep["positives"]), ", ".join(sep["nematode_families"]), len(sep["negatives"]), len({(R[i].get("extra") or {}).get("family") for i in sep["negatives"]}), cs["candidates"], {k: round(v, 2) for k, v in cs["auc_by_family"].items()}, cs["auc_mean"], cs["n_null"], cs["null_mean"], cs["null_p95"], cs["p"], ct["candidates"], ct["auc_mean"], ct["null_mean"], ct["null_p95"], ct["p"]))
open(os.path.join(CRY, "docs", "domain_tree_discordance.md"), "w", encoding="utf-8").write("""# Domain trees compared entry by entry (2026-11-27)

For every entry the 10 nearest relatives (patristic distance) in the Domain I, II and III trees (structure-based alignments); the majority group of their names at the secondary rank (Cry1A, Cry1B ...) is the yardstick.
Entries whose nearest relatives are not dominated (>= 60 %%) by one group in some tree cannot be judged: %d of %d (mostly the small divergent families whose names do not form groups).
Of the %d assessable entries: concordant %d, one domain differs %d (Domain I %d, II %d, III %d), all three differ %d.
Noise: the same domain from the structure-based and the sequence-based alignment gives a different majority group in %s of the entries with a clear neighbourhood; between different domains the rate is 2.6-3.7 %%.
Normalised Robinson-Foulds distance between the domain trees: %s; between the two alignments of the same domain (noise): %s. (RF is sensitive to the many unresolved near-identical sequences.)

Entries with a domain that differs (Domain III dominates, in Cry1B, Cry1C, Cry1E and Cry1J: the Domain III exchange known from the Cry1 family):
%s
Full table: domain_tree_discordance.tsv. Script 25_domain_trees.py.
""" % (cc.get("mixed neighbourhood", 0), dd["n_entries"], assess, cc.get("concordant", 0), cc.get("Domain I differs", 0) + cc.get("Domain II differs", 0) + cc.get("Domain III differs", 0), cc.get("Domain I differs", 0), cc.get("Domain II differs", 0), cc.get("Domain III differs", 0), cc.get("all three differ", 0),
       {k: "%d of %d" % tuple(v) for k, v in dd["replicate_noise"].items()}, {k: round(v, 2) for k, v in dd["rf"].items()}, {k: round(v, 2) for k, v in dd["rf_replicate"].items()},
       "\n".join("* %s: %s (groups DI / DII / DIII = %s)" % (R[x]["label"].split(" ")[0], dd["category"][x], " / ".join(dd["groups"][x][d][0] for d in ("DI", "DII", "DIII"))) for x in top)))
print("activity:", dict(cnt_a)); print("domain trees:", dict(cnt_d2)); print("marks:", [(m["id"], len(m["cols"])) for m in meta["marks"]]); print("tangle", meta["tangle"])
