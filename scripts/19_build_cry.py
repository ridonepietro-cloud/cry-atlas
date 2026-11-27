"""Build the Cry pack (site/data/cry/): rows (work/rows.json), frames, modules, classes, specificity text -> cry_pack.build().
Frames: 'struct' = FoldMason structure-based frame (work/aln/struct.fa, when it exists) and 'seq' = segment-wise MAFFT frame (work/aln/seqseg.fa).  The first frame is the default (id 'struct' in meta).
Modules = N-term | DI | DII | DIII | C-term: the segments of the MAFFT frame; for the structure frame the column ranges come from the module of every residue (monotone dynamic programme over the columns).
Run: py -I 19_build_cry.py [--no-struct]"""
import os, sys, re, json, collections, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
import cry_pack as CP
rows = jload(os.path.join(WORK, "rows.json")); classes = jload(os.path.join(WORK, "classes.json"), {})
seg = CC.read_fasta(os.path.join(WORK, "aln", "seqseg.fa")); seg = dict(seg)
segm = jload(os.path.join(WORK, "aln", "seqseg_modules.json"))["modules"]; NAMES = [m[0] for m in segm]
use_struct = os.path.exists(os.path.join(WORK, "aln", "structseg.fa")) and "--no-struct" not in sys.argv
# per-residue module (from the segment frame)
col2mod = {}
for n, a, b in segm:
    for c in range(a, b + 1): col2mod[c] = NAMES.index(n)
resmod = {r["id"]: [col2mod[c] for c in CC.col_index(seg[r["id"]])] for r in rows}


def modules_for(gapped):
    """column ranges of the modules in another frame: monotone DP over the columns, votes = residues of the rows with that module at that column"""
    W = len(next(iter(gapped.values()))); votes = [[0] * len(NAMES) for _ in range(W + 1)]
    for rid, g in gapped.items():
        for c, m in zip(CC.col_index(g), resmod[rid]): votes[c][m] += 1
    NEG = -10 ** 9; K = len(NAMES); best = [[NEG] * K for _ in range(W + 1)]; back = [[0] * K for _ in range(W + 1)]
    for s in range(K): best[0][s] = 0
    for c in range(1, W + 1):
        run = (NEG, 0)
        for s in range(K):
            run = max(run, (best[c - 1][s], s)); best[c][s] = run[0] + votes[c][s]; back[c][s] = run[1]
    s = max(range(K), key=lambda t: best[W][t]); path = [0] * (W + 1)
    for c in range(W, 0, -1): path[c] = s; s = back[c][s]
    mods = []
    for c in range(1, W + 1):
        if mods and mods[-1][0] == path[c]: mods[-1][2] = c
        else: mods.append([path[c], c, c])
    return [(NAMES[k], a, b) for k, a, b in mods]


frames = []
if use_struct:
    st = dict(CC.read_fasta(os.path.join(WORK, "aln", "structseg.fa")))
    frames.append({"id": "struct", "label": "Structure-based: FoldMason of the AlphaFold DB models of each module (N-term, DI, DII, DIII, C-term cut separately), the entries without an exact model added with MAFFT --add, segments joined", "kind": "structure", "short": "Structure-based", "gapped": st})
    modules = [tuple(m) for m in jload(os.path.join(WORK, "aln", "structseg_modules.json"))["modules"]]
    frames.append({"id": "seq", "label": "Sequence-based: MAFFT L-INS-i of N-term, DI, DII, DIII and C-term separately, joined (segments from the first-pass alignment)", "kind": "sequence", "short": "Sequence-based", "gapped": seg})
else:
    frames.append({"id": "struct", "label": "Sequence-based: MAFFT L-INS-i of N-term, DI, DII, DIII and C-term separately, joined (segments from the first-pass alignment)", "kind": "sequence", "short": "Sequence-based", "gapped": seg})
    modules = [tuple(m) for m in segm]
print("frames", [(f["id"], len(next(iter(f["gapped"].values())))) for f in frames], "modules", modules, flush=True)
MODCOL = ["N-term", "DI", "DII", "DIII", "C-term"]


def short_species(sp):
    p = sp.split()
    if not p: return ""
    if p[0] == "Bacillus" and len(p) > 1: return "B. " + " ".join(p[1:3] if p[1] != "thuringiensis" else ["thuringiensis"] + (["sv. " + p[3]] if len(p) > 3 and p[2] == "serovar" else []))
    return " ".join(p[:2])


def tgt_text(t):
    if not t: return None
    bits = []
    if t["yes"]: bits.append("active against " + ", ".join("%s (%d)" % (o, c) for o, c in sorted(t["yes"].items(), key=lambda x: -x[1])) + ": " + ", ".join(t["sp_yes"][:8]) + (" ..." if len(t["sp_yes"]) > 8 else ""))
    if t["no"]: bits.append("no activity recorded against " + ", ".join("%s (%d)" % (o, c) for o, c in sorted(t["no"].items(), key=lambda x: -x[1])) + ": " + ", ".join(t["sp_no"][:6]) + (" ..." if len(t["sp_no"]) > 6 else ""))
    return "BPPRC specificity data: " + "; ".join(bits)


for r in rows:
    r["class"] = classes.get(r["family"], "other")
    r["gene"] = r["name"].replace("Cry", "cry", 1)
    r["acc"] = r["accessions"][0][1]
    sp = short_species(r["species"]); more = (" +%d" % (len(r["names"]) - 1)) if len(r["names"]) > 1 else ""
    r["label"] = "%s%s · %s" % (r["name"], more, sp) if sp else "%s%s" % (r["name"], more)
    r["header"] = "%s | %s | %s" % (r["name"], r["acc"], r["species"] or "unknown source")
    info = []
    info.append("BPPRC Cry entry%s: %s" % ("" if len(r["names"]) == 1 else " (%d names, identical sequence)" % len(r["names"]), ", ".join(r["names"])))
    if r["old_names"] or r["other_names"]: info.append("Other names: " + ", ".join(r["old_names"] + r["other_names"]))
    info.append("BPPRC year %s%s" % (r["year"] or "n/a", "; holotype sequence of its name in BPPRC" if r["holotype"] else ""))
    d = r["bpprc_domains"]
    if d.get("d1") or d.get("d2") or d.get("d3"): info.append("BPPRC domain boundaries (residues of this sequence): " + "; ".join("D%d %s" % (k, d.get("d%d" % k)) for k in (1, 2, 3) if d.get("d%d" % k)))
    else: info.append("BPPRC gives no domain boundaries for this entry (modules here come from the alignment)")
    if r["strains"] or r["plasmids"] or r["species"]: info.append("Source: %s%s%s" % (r["species"] or "not stated", "; strain " + ", ".join(r["strains"][:3]) if r["strains"] else "", "; plasmid " + ", ".join(r["plasmids"][:2]) if r["plasmids"] else ""))
    if r.get("uniprot"): info.append("UniProt: " + ", ".join(r["uniprot"][:4]))
    tt = tgt_text(r.get("target"))
    if tt: info.append(tt)
    r["info"] = info + r["info"]
    t = r.get("target") or {}
    r["extra"] = {"names": r["names"], "family": r["family"], "year": r["year"], "holotype": r["holotype"], "target": sorted(t.get("yes", {}), key=lambda o: -t["yes"][o]), "nontarget": sorted(t.get("no", {}))}

# ---- colour-by attributes (page: colour entries by ...; filter chips) ----
PAL = ["#d64045", "#3b7cbd", "#2f8a57", "#d49a14", "#8a5cc4", "#1f9d8f", "#c2495f", "#8a94a0", "#b0793a"]


def serovar(sp):
    m = re.search(r"serovar (\w+)", sp or "")
    return m.group(1) if m else None


def a_target(r):
    t = r.get("target")
    if t and t["yes"]:
        o = max(t["yes"].items(), key=lambda x: x[1])[0]; return "Nematoda (Rhabditida)" if o == "Rhabditida" else o if o in ("Lepidoptera", "Coleoptera", "Diptera", "Hemiptera") else "other order"
    return "none active recorded" if t else "no data"


def a_era(r):
    y = r["year"]; return "n/a" if not y else "to 1995" if y <= 1995 else "1996-2005" if y <= 2005 else "2006-2015" if y <= 2015 else "2016 or later"


def a_source(r):
    sp = r["species"] or ""
    if not sp or sp.startswith("Unknown"): return "unknown"
    if "thuringiensis" in sp:
        sv = serovar(sp); return "Bt " + sv if sv in ("kurstaki", "aizawai", "israelensis", "entomocidus") else "Bt other serovar" if sv else "Bt (no serovar stated)"
    return "other species"


def a_model(r):
    ms = r["models"]
    if not ms: return "no model"
    if any("donor" in m for m in ms): return "model of a near-identical entry"
    exp = any(m["kind"] == "experimental" for m in ms); af = any(m["kind"] == "AF2" for m in ms)
    return "PDB + AlphaFold DB" if exp and af else "PDB only" if exp else "AlphaFold DB"


def a_len(r):
    return "long (>= 800 aa, protoxin)" if r["len"] >= 800 else "short (< 800 aa)"


ATTRS = [("target", "Target order (BPPRC specificity)", "Target order", "Insect order (or nematodes) with most BPPRC activity records; 'none active recorded' = only negative records; 'no data' = no specificity record for any of the names", a_target),
         ("host", "Source organism (GenPept)", "Source", "Organism of the GenPept record; Bt serovars named in the record", a_source),
         ("era", "Year of the BPPRC entry", "Year", "Year given by BPPRC for the first name of the sequence", a_era),
         ("size", "Length (protoxin or toxin only)", "Length", "Sequences of 800 residues or more carry the C-terminal protoxin half", a_len),
         ("model", "Structure model", "Model", "Source of the models of the entry (primary model = most residues placed)", a_model)]
attrs_def = []
for key, label, short, text, fn in ATTRS:
    cnt = collections.Counter(fn(r) for r in rows)
    vals = [v for v, _ in sorted(cnt.items(), key=lambda x: (-x[1]))]
    for r in rows: r.setdefault("attrs", {})[key] = fn(r)
    vcol = {v: PAL[i] if i < len(PAL) - 1 else PAL[-1] for i, v in enumerate(vals)}
    for v in ("no data", "n/a", "unknown", "no model", "none active recorded"):
        if v in vcol: vcol[v] = "#b9bec4"
    attrs_def.append({"key": key, "label": label, "short": short, "text": text, "values": [{"id": v, "label": v, "colour": vcol[v], "n": cnt[v]} for v in vals]})
features = {"attrs": attrs_def, "arch": [], "arch_sub": [], "sgroups": [], "note": "attributes of the Cry rows for the 'Colour entries by' bar"}

ids = [r["id"] for r in rows]
order = jload(os.path.join(WORK, "aln", "seqseg_order.json"))["order"]
ref = next(r["id"] for r in rows if r["name"] == "Cry1Aa1")
spec = {"pack": "cry", "system": "cry", "title": "Cry toxins (3-domain insecticidal proteins of the BPPRC database)",
        "build": "Cry pack built %s from the BPPRC Cry class (bpprc-db.org, 797 sequences downloaded 2026-11-24; 670 distinct sequences); AlphaFold DB models + PDB entries; RING 4.0" % time.strftime("%Y-%m-%d"),
        "seg": "Cry toxin", "frames": frames, "order": order, "rows": rows, "modules": modules, "ref_row": ref,
        "default": {"auto": True, "row": ref, "why": "inter-domain NODE column of Cry1Aa1 seen in the most models"},
        "node_text": "NODE (adapted rule of the packages): a column is a strict NODE when the contacts that join two modules (domains) at that column, over all models of the pack, have >= 2 interaction types, one not VDW, and are seen in >= %d models; "
                     "a residue is NODE when it makes such an inter-domain contact in the row's own primary model at a strict column; the pale mark = residue at a strict column in a row that does not show the contact (no model, or no contact there)." % CP.NODE_MIN,
        "node_short": "yes (contact with another domain at a strict NODE column)",
        "order_basis": "guide-tree order of the Domain II alignment (MAFFT --reorder); the ML tree order is used when the Trees view is on", "order_short": "tree",
        "inputs": {"bpprc": "bpprc-db.org Cry class (category_download 2026-11-24), holotype list, per-domain sequences, structures table, specificity export", "structures": "AlphaFold DB (exact-sequence models), RCSB PDB entries >= 97 % identical to a row", "contacts": "RING 4.0"},
        "cry": {"modules": MODCOL, "n_names": sum(len(r["names"]) for r in rows), "features": {"file": "features.json"}}, "features": features, "tangle": {"default": ["struct_DII", "struct_DIII"]}}
CP.build(spec)
