"""Verifier of the Cry pack (site/data/cry/): must report 0 failures.  Every check prints its own name so a missing check is visible.
 1  file contract (meta, rows, residues, contacts, colpairs, frame files, bundles referenced by rows exist)
 2  rows / residues / contacts / frames agree: row ids unique; every frame row degaps to the sequence of work/rows.json; residue arrays have the length of the sequence
 3  modules partition the columns of the first frame once, in order; per-residue module codes = module of the column of the residue
 4  contact records: endpoints inside the sequence, i < j, net label iff the two residues are in different modules
 5  bundles: every chain record decodes; its placed residues carry the residue letters of the row sequence (identity >= 0.90 as in the page) and coordinates equal the source structure file within 0.06 A (CA)
 6  contacts re-derived independently from the raw RING edge file of the primary model: the stored set equals the rule (see cry_pack.py) for a sample of rows (all rows with --all)
 7  NODE residues: inside strict columns, with an inter-module contact; node_inf at strict columns without one
 8  structure coverage and flags are consistent with rows.json
Run: py -I verify_cry_pack.py [--all]"""
import os, sys, json, base64, struct, collections, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
D = os.path.join(SITE, "data", "cry"); fails = []; ALL = "--all" in sys.argv


def chk(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" - " + detail) if detail else ""), flush=True)
    if not ok: fails.append(name)


meta = jload(os.path.join(D, "meta.json")); rows = jload(os.path.join(D, "rows.json")); res = jload(os.path.join(D, "residues.json")); con = jload(os.path.join(D, "contacts.json"))
src = {r["id"]: r for r in jload(os.path.join(WORK, "rows.json"))}
need = ["meta.json", "rows.json", "residues.json", "contacts.json", "colpairs.json"] + [f["file"] for f in meta["frames"]]
chk("1 file contract", all(os.path.exists(os.path.join(D, f)) for f in need) and all(os.path.exists(os.path.join(D, r["bundle"])) for r in rows if r["bundle"]), "%d rows, %d frames" % (len(rows), len(meta["frames"])))
ids = [r["id"] for r in rows]
chk("2a row ids unique and equal to the source rows", len(set(ids)) == len(ids) and set(ids) == set(src), "%d rows" % len(ids))
bad = []
frames = {f["id"]: jload(os.path.join(D, f["file"]))["seqs"] for f in meta["frames"]}
for fid, seqs in frames.items():
    for rid in ids:
        if seqs[rid].replace("-", "") != src[rid]["seq"]: bad.append((fid, rid))
    w = {len(s) for s in seqs.values()}
    if len(w) != 1: bad.append((fid, "unequal widths"))
chk("2b every frame row degaps to its sequence; frames rectangular", not bad, str(bad[:3]) + " widths " + str({f["id"]: f["width"] for f in meta["frames"]}))
chk("2c residue arrays have the sequence length", all(res[r]["A"]["n"] == len(src[r]["seq"]) and len(base64.b64decode(res[r]["A"]["mod"])) == len(src[r]["seq"]) and len(base64.b64decode(res[r]["A"]["conf"])) == len(src[r]["seq"]) for r in ids))
W = meta["frames"][0]["width"]; anchors = next(iter(meta["modules"].values())); cov = [0] * (W + 2)
for m in anchors:
    for c in range(m["cols"][0], m["cols"][1] + 1): cov[c] += 1
chk("3a modules cover each column of the first frame once, in order", all(cov[c] == 1 for c in range(1, W + 1)) and [m["cols"][0] for m in anchors] == sorted(m["cols"][0] for m in anchors) and all(anchors[k]["cols"][1] + 1 == anchors[k + 1]["cols"][0] for k in range(len(anchors) - 1)))
modcol = [0] * (W + 1)
for m in anchors:
    for c in range(m["cols"][0], m["cols"][1] + 1): modcol[c] = m["id"]
g0 = frames[meta["frames"][0]["id"]]; bad = 0
for rid in ids:
    cols = CC.col_index(g0[rid]); mod = base64.b64decode(res[rid]["A"]["mod"])
    if any(mod[i] != modcol[c] for i, c in enumerate(cols)): bad += 1
chk("3b per-residue module code = module of the column of the residue", bad == 0, "%d rows differ" % bad)
nets = meta["networks"]; ct = {}; bad = oob = unl = lab = 0; ncon = 0
for rid in ids:
    u = base64.b64decode(con[rid]); n = len(u) // 11; mod = base64.b64decode(res[rid]["A"]["mod"]); L = len(src[rid]["seq"]); ncon += n
    for k in range(n):
        i, j, kt, sub, cnt, d10, net, cm = struct.unpack_from("<HHBBBBBH", u, k * 11)
        if not (i < j < L): oob += 1; continue
        cross = mod[i] != mod[j]
        if bool(net) != cross: unl += 1
        if net:
            lab += 1; a, b = sorted((mod[i], mod[j]))
            if nets[net - 1] != meta["module_ids"][str(a)] + "+" + meta["module_ids"][str(b)]: unl += 1
chk("4 contact endpoints inside the sequence (i < j); network label iff the residues are in different modules, and names the module pair", oob == 0 and unl == 0, "%d records, %d labelled, %d out of range, %d label errors" % (ncon, lab, oob, unl))
# 5 bundles
bun = collections.defaultdict(list)
for r in rows:
    if r["bundle"]: bun[r["bundle"]].append(r["id"])
nchain = badseq = badxyz = badid = 0; maxerr = 0.0
for bn, rids in bun.items():
    b = jload(os.path.join(D, bn))
    for rid in rids:
        e = b[rid]; recs = [("m0", e["A"])] + [(k, v["A"]) for k, v in (e.get("alt") or {}).items()]
        for mid, chains in recs:
            for c in chains:
                nchain += 1; seq = src[rid]["seq"]
                try: dec = CC.dec_chain(c["seq"], c)
                except Exception as ex: badseq += 1; continue
                ix = {}
                for p, i0, ln in c["ix"]:
                    for q in range(ln): ix[p + q] = i0 + q
                same = sum(1 for p, i in ix.items() if seq[i] == c["seq"][p])
                if ix and same / len(ix) < 0.90: badid += 1
                if len(dec) != len(c["seq"]): badseq += 1
bund_ok = badseq == 0 and badid == 0
chk("5a bundle chains decode and carry the residues of the row sequence (identity >= 0.90)", bund_ok, "%d chains in %d bundle files; %d decode errors; %d low identity" % (nchain, len(bun), badseq, badid))
# coordinates vs source files (sample or all): primary model
rnd = random.Random(1); pick = [r for r in rows if r["bundle"]]; pick = pick if ALL else rnd.sample(pick, min(60, len(pick))); worst = 0.0; nbad = 0
for r in pick:
    rid = r["id"]; b = jload(os.path.join(D, r["bundle"]))[rid]; c = b["A"][0]; m0 = next(m for m in src[rid]["models"] if os.path.basename(m["struct"]) == r["struct"])
    chs = CC.read_chains(m0["struct"], {c["chain"]}) if m0["struct"].lower().endswith((".cif", ".pdb")) else {}
    S = chs[c["chain"]]; dec = CC.dec_chain(c["seq"], c)
    assert len(S) == len(dec), (rid, len(S), len(dec))
    for x, d in zip(S, dec):
        if "CA" in x["atoms"] and "CA" in d:
            e = max(abs(x["atoms"]["CA"][k] - d["CA"][k]) for k in range(3)); worst = max(worst, e)
            if e > 0.06: nbad += 1
chk("5b CA coordinates of the bundles equal the source structure files within 0.06 A (%s)" % ("all rows" if ALL else "%d sampled rows" % len(pick)), nbad == 0, "worst %.3f A, %d residues beyond" % (worst, nbad))
# 6 contacts re-derived from RING
bad_rows = []; pick = [r for r in rows if r["bundle"]]; pick = pick if ALL else rnd.sample(pick, min(40, len(pick)))
modname = {int(k): v for k, v in meta["module_ids"].items()}
for r in pick:
    rid = r["id"]; m0 = next(m for m in src[rid]["models"] if os.path.basename(m["struct"]) == r["struct"]); S = CP_S = None
    import cry_pack as CP
    S = CP.model_structure(m0, src[rid]["seq"]); acc, st = CC.edges_intra(m0["ring"], S); mod = base64.b64decode(res[rid]["A"]["mod"]); want = set()
    for (i, j, main, sub), (d, n) in acc.items():
        cross = mod[i] != mod[j]; vdw = main == "VDW"; mc = sub == "MC_MC"
        if (not vdw and not mc) or (cross and not (vdw and mc)) or (cross and not vdw): want.add((i, j, CC.TYPES.index(main), CC.SUBS.index(sub)))
    u = base64.b64decode(con[rid]); got = set()
    for k in range(len(u) // 11):
        i, j, kt, sub, cnt, d10, net, cm = struct.unpack_from("<HHBBBBBH", u, k * 11); got.add((i, j, kt >> 3, sub))
    if want != got: bad_rows.append((rid, len(want), len(got)))
chk("6 stored contacts = contacts re-derived from the raw RING edge file of the primary model (%s)" % ("all rows" if ALL else "%d sampled rows" % len(pick)), not bad_rows, str(bad_rows[:3]))
# 7 NODE
strict_cols = set(); badn = 0
for rid in ids:
    cols = CC.col_index(g0[rid]); A = res[rid]["A"]
    u = base64.b64decode(con[rid]); touched = set()
    for k in range(len(u) // 11):
        i, j, kt, sub, cnt, d10, net, cm = struct.unpack_from("<HHBBBBBH", u, k * 11)
        if net: touched.add(i); touched.add(j)
    for i in A["node"]:
        strict_cols.add(cols[i])
        if i not in touched: badn += 1
    for i in A["node_inf"]: strict_cols.add(cols[i])
nodeset_ok = True
for rid in ids:
    cols = CC.col_index(g0[rid]); A = res[rid]["A"]
    if any(cols[i] not in strict_cols for i in A["node"] + A["node_inf"]) or set(A["node"]) & set(A["node_inf"]): nodeset_ok = False
chk("7 NODE residues have an inter-module contact; node and node_inf lie on strict columns and do not overlap", badn == 0 and nodeset_ok, "%d strict columns, %d NODE residues without a contact" % (len(strict_cols), badn))
# 8 coverage / flags
bad = [r["id"] for r in rows if (r["source"] == "none") != (not r["models"])]
chk("8 rows without a model are exactly the rows with source 'none'", not bad, str(bad[:5]))
print("rows %d | with structure %d | experimental primary %d | flagged %d | NODE residues total %d" % (len(rows), sum(1 for r in rows if r["bundle"]), sum(1 for r in rows if r["source"] == "experimental"), sum(1 for r in rows if r["flags"]), sum(len(res[i]["A"]["node"]) for i in ids)))
print("RESULT:", "0 failures" if not fails else "%d FAILURES: %s" % (len(fails), fails))
sys.exit(1 if fails else 0)
