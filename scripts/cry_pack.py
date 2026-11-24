"""Pack writer of the Cry Atlas (adaptation of homo_pack.py of the NODE-MSA Atlas for a MONOMERIC multi-domain protein): rows, residues, contacts, structure bundles, inferred layers  ->  site/data/<pack>/

File contract of the Atlas packs (meta, rows, aln_*, residues, contacts, bundle_NNN, colpairs); subunit "A" is the toxin chain, there is no subunit B (meta.homo = true).
What differs from the oligomer packs: the toxins are single-chain proteins, so the 'interface' of the NODE rule is the contact between MODULES (domains) inside the chain:
  * contact kind 0 (same chain) for every stored contact; net = pair of modules joined by the contact (D1+D2, D2+D3, ...), 0 for a contact inside one module;
  * strict NODE column (rule of the packages at column level, applied to inter-module contacts): over all models of the pack the contacts that join two modules at that column have >= 2 interaction types,
    one of them not VDW, and are seen in >= NODE_MIN models; a residue is NODE when it makes an inter-module contact in the row's own primary model at such a column; node_inf = residue at such a column
    in a row that does not show it (no model, or no inter-module contact there);
  * stored contacts (size limit of the page): all side-chain hydrogen bonds and other non-VDW, non-backbone contacts, and every inter-module contact except backbone-backbone VDW; VDW contacts inside a module are not stored.
spec = { pack, system, title, build, seg, frames [ {id, label, kind, short, gapped {row id: str}} ] (first = default frame 'struct'), order [row ids], rows [...], modules [(name, first col, last col)] in the columns of the first frame,
         ref_row, default {...}, inputs {...}, node_text }
row dict: id, label, species, gene, class, acc, accs, pdb, header, seq, flags [], info [], models [ {label, kind, struct, chains, ring, pdb, ...} ]
"""
import os, sys, json, collections, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import SITE
import cry_core as CC
from cry_core import TYPES, SUBS, col_index, b64, pack_contacts, chain_records, runs_of

BUNDLE_ROWS = 6
NODE_MIN = 5                                                     # a strict NODE column needs inter-module contacts in at least this many models


def model_structure(m, seq):
    """chain dict of the protomer of a model, placed on the row sequence; the chain is the first of m['chains'] present in the file with >= 50 residues"""
    path = m["struct"]
    chs = CC.read_chains(path, set(m["chains"]) if m["kind"] == "experimental" else None)
    ch = next((c for c in m["chains"] if c in chs and len(chs[c]) >= 50), None) if m["kind"] == "experimental" else next(iter(chs))
    if ch is None:
        raise ValueError("%s: none of the chains %s with >= 50 residues" % (path, m["chains"]))
    scheme = CC.cif_scheme(path) if path.lower().endswith((".cif", ".mmcif")) else None
    return CC.one_chain_structure(path, ch, seq, scheme)


def build(spec, verbose=True):
    t0 = time.time()
    PACK = spec["pack"]; OUT = os.path.join(SITE, "data", PACK)
    rows = spec["rows"]; frames = spec["frames"]; F0 = frames[0]
    W = len(next(iter(F0["gapped"].values())))
    MODULES = spec["modules"]; MOD_ID = {m[0]: k + 1 for k, m in enumerate(MODULES)}; MOD_NAME = {v: k for k, v in MOD_ID.items()}
    assert len(MODULES) <= 10, "the page has 10 module colours"
    modcol = [0] * (W + 1)
    for name, a, b in MODULES:
        for c in range(a, b + 1):
            modcol[c] = MOD_ID[name]
    NETS = {}
    for r in rows:
        for F in frames:
            assert F["gapped"][r["id"]].replace("-", "") == r["seq"], "row %s: frame %s differs from the sequence" % (r["id"], F["id"])
        r["cols"] = col_index(F0["gapped"][r["id"]])

    classes_all = sorted({r["class"] for r in rows})
    ctypes = collections.defaultdict(set); csid = collections.defaultdict(set); nmodels_at = collections.Counter()
    built = {}; cpairs = collections.defaultdict(lambda: collections.defaultdict(set)); nstruct = collections.Counter(); stats = collections.Counter()

    def model_records(r, m):
        """(S, stored records, interface info) of one model of a row; contacts that join two modules are labelled with the module pair"""
        S = model_structure(m, r["seq"]); cols = r["cols"]
        acc, st = CC.edges_intra(m["ring"], S) if m.get("ring") and os.path.exists(m["ring"]) else ({}, collections.Counter())
        recs, inter = [], []
        for (i, j, main, sub), (d, n) in acc.items():
            ma, mb = modcol[cols[i]], modcol[cols[j]]
            cross = ma != mb
            if cross:
                a, b = sorted((ma, mb)); net = NETS.setdefault(MOD_NAME[a] + "+" + MOD_NAME[b], len(NETS) + 1)
                inter.append((i, j, main, sub))
            else:
                net = 0
            vdw = main == "VDW"; mc = sub == "MC_MC"
            keep = (not vdw and not mc) or (cross and not (vdw and mc)) or (cross and not vdw)
            if keep:
                recs.append([i, j, 0, TYPES.index(main), SUBS.index(sub), min(n, 255), int(round(d * 10)), net, 1])
        return S, recs, inter, st

    for r in rows:
        flags, info = list(r.get("flags", [])), list(r.get("info", []))
        outs = []
        for m in r["models"]:
            try:
                S, recs, inter, st = model_records(r, m)
            except Exception as e:                                                           # a model that cannot be read is reported, not hidden
                flags.append("model %s could not be used: %s" % (m.get("label", "?"), str(e)[:120])); stats["model_error"] += 1; continue
            c = S["chains"][S["rep_A"]]
            if c["ident"] < 0.9 or c["placed"] < 0.5 * len(c["res"]):
                flags.append("%s: the model residues agree with the sequence at %.0f %% (%d of %d residues placed): not used" % (m["label"], 100 * c["ident"], c["placed"], len(c["res"]))); stats["model_mismatch"] += 1; continue
            if m["kind"] == "experimental" and c["ident"] < 0.999:
                info.append("%s differs from the sequence of this entry at %d placed residue(s) (mutant or other variant: %.1f %% identity)" % (m["pdb"], round((1 - c["ident"]) * c["placed"]), 100 * c["ident"]))
            outs.append({"m": m, "S": S, "recs": recs, "inter": inter, "placed": c["placed"], "ident": c["ident"]})
        # primary = the model that places most residues of the row sequence (experimental wins a tie): full-length first
        outs.sort(key=lambda o: (-o["placed"], 0 if o["m"]["kind"] == "experimental" else 1))
        for k, o in enumerate(outs):
            mid = r["id"] + "|m%d" % k; nstruct[r["class"]] += 1
            for i, j, main, sub in o["inter"]:
                for c in (r["cols"][i], r["cols"][j]):
                    ctypes[c].add(main); csid[c].add(mid)
            for i, j, kd, typ, sb, n, d10, net, cm in o["recs"]:
                if net:
                    ca, cb = r["cols"][i], r["cols"][j]; cpairs[r["class"]][(min(ca, cb), max(ca, cb), 1)].add(mid)
        if outs and outs[0]["placed"] < len(r["seq"]):
            info.append("the primary model covers %d of %d residues of the sequence" % (outs[0]["placed"], len(r["seq"])))
        if not outs:
            if not any("could not be used" in f or "not used" in f for f in flags):
                info.append("No structure model for this entry (sequence, modules and the inferred layers only).")
        built[r["id"]] = (outs, flags, info)
    strict = {c for c, ts in ctypes.items() if len(ts) >= 2 and any(t != "VDW" for t in ts) and len(csid[c]) >= NODE_MIN}
    if verbose:
        print("strict NODE columns:", len(strict), "of", W, "columns | models used", sum(nstruct.values()), "| stats", dict(stats), "| networks", list(NETS), flush=True)

    res_out, contacts_out, bundles, rows_out = {}, {}, {}, []
    for r in rows:
        n = len(r["seq"]); cols = r["cols"]; outs, flags, info = built[r["id"]]
        mod = [modcol[c] for c in cols]; pl = [255] * n
        o0 = outs[0] if outs else None
        if o0:
            c = o0["S"]["chains"][o0["S"]["rep_A"]]
            for x, i in zip(c["res"], c["idx"]):
                if i >= 0:
                    pl[i] = 254 if o0["m"]["kind"] == "experimental" else max(0, min(100, int(round(x["bfac"]))))
        touched = set()
        if o0:
            for q in o0["recs"]:
                if q[7]:
                    touched.add(q[0]); touched.add(q[1])
        nodes = sorted(i for i in touched if cols[i] in strict); ns = set(nodes)
        inf = [i for i in range(n) if i not in ns and cols[i] in strict]
        res = {"A": {"n": n, "mod": b64(mod), "conf": b64(pl), "node": nodes, "node_inf": inf, "chains": 1 if o0 else 0},
               "B": {"n": 0, "mod": "", "conf": "", "node": [], "node_inf": [], "chains": 0}, "pkg_node": {"A": [], "B": []}}
        bund = {}; models = []; alt_out = {}
        if o0:
            S0 = o0["S"]; bund["A"] = chain_records(S0, S0["order"]["A"]); bund["B"] = []; c = S0["chains"][S0["rep_A"]]
            res["A"].update({"rep": S0["rep_A"], "chain_len": len(c["res"]), "placed": c["placed"], "ident": round(c["ident"], 4)})
            res["A"]["rn"] = runs_of(sorted(((i, x["resnum"]) for x, i in zip(c["res"], c["idx"]) if i >= 0)))
            models.append({"id": "m0", "label": o0["m"]["label"], "kind": o0["m"]["kind"], "file": os.path.basename(o0["m"]["struct"]), "heptamer": False, "primary": True})
        for k, o in enumerate(outs[1:], 1):
            mid = "m%d" % k
            alt_out[mid] = {"A": chain_records(o["S"], o["S"]["order"]["A"]), "B": [], "con": pack_contacts(o["recs"]), "chainsA": o["S"]["order"]["A"], "chainsB": []}
            models.append({"id": mid, "label": o["m"]["label"], "kind": o["m"]["kind"], "file": os.path.basename(o["m"]["struct"]), "heptamer": False, "primary": False, "n_contacts": len(o["recs"])})
        if alt_out:
            bund["alt"] = alt_out
        res_out[r["id"]] = res; contacts_out[r["id"]] = pack_contacts(o0["recs"] if o0 else []); bundles[r["id"]] = bund
        pe = next((o["m"] for o in outs if o["m"]["kind"] == "experimental"), None)
        rows_out.append({"id": r["id"], "label": r["label"], "species": r["species"], "gene": r["gene"], "acc": r["acc"], "accB": "", "pdb": pe["pdb"] if pe else None, "header": r["header"], "class": r["class"],
                         "pair": r["id"], "source": ("experimental" if o0 and o0["m"]["kind"] == "experimental" else "AF2 monomers" if o0 else "none"),
                         "struct": os.path.basename(o0["m"]["struct"]) if o0 else "", "chainsA": o0["S"]["order"]["A"] if o0 else [], "chainsB": [], "nA": n, "nB": 0, "flags": flags, "info": info,
                         "n_contacts": len(o0["recs"]) if o0 else 0, "n_labelled": sum(1 for q in o0["recs"] if q[7]) if o0 else 0, "models": models,
                         "extra": r.get("extra", {}), **r.get("attrs", {})})
    colpairs = {c: sorted([a, b, len(v), cl] for (a, b, cl), v in d.items() if len(v) >= 3) for c, d in cpairs.items()}
    ref = next(x for x in rows if x["id"] == spec["ref_row"]); anchors = []
    for name, a, b in MODULES:
        inside = [i for i, c in enumerate(ref["cols"]) if a <= c <= b]
        anchors.append({"id": MOD_ID[name], "name": name, "cols": [a, b], "ref_row": ref["id"], "ref_sub": "A", "ref_idx": [inside[0], inside[-1]] if inside else None})
    dflt = spec["default"]
    if dflt.get("auto"):                                                     # the NODE column seen in inter-module contacts of most models, at the reference row
        have = set(ref["cols"]); best = max(((len(csid[c]), c) for c in strict if c in have), default=None)
        dflt = {"row": dflt["row"], "col": best[1] if best else ref["cols"][len(ref["cols"]) // 2], "why": dflt["why"] + (" (inter-domain contacts in %d models)" % best[0] if best else "")}
    os.makedirs(OUT, exist_ok=True)
    for f in os.listdir(OUT):
        if f in ("gc.json", "gcsyn.json", "trees.json", "specificity.json") or f.startswith("asr_") or f.startswith("species"):
            continue
        if f.startswith("bundle_") or f.startswith("aln_") or f.endswith(".json"):
            os.remove(os.path.join(OUT, f))

    def dump(name, obj):
        p = os.path.join(OUT, name); json.dump(obj, open(p, "w", encoding="utf-8"), separators=(",", ":")); return os.path.getsize(p)
    sizes = {}
    order = spec["order"]; assert sorted(order) == sorted(x["id"] for x in rows_out), "order list differs from the rows"
    rid_order = {i: k for k, i in enumerate(order)}; rows_out.sort(key=lambda x: rid_order[x["id"]]); ids_b = [x["id"] for x in rows_out]
    for k, g in enumerate([ids_b[i:i + BUNDLE_ROWS] for i in range(0, len(ids_b), BUNDLE_ROWS)]):
        name = "bundle_%03d.json" % k
        if not any(bundles[i] for i in g):
            continue
        sizes[name] = dump(name, {i: bundles[i] for i in g})
        for ro in rows_out:
            if ro["id"] in g:
                ro["bundle"] = name if bundles[ro["id"]] else None
    for ro in rows_out:
        ro.setdefault("bundle", None)
    classes = sorted({r["class"] for r in rows}); minsup = {c: 3 for c in classes}
    fr_meta = []
    for k, F in enumerate(frames):
        fname = "aln_structural.json" if k == 0 else "aln_%s.json" % F["id"]
        sizes[fname] = dump(fname, {"seqs": {r["id"]: F["gapped"][r["id"]] for r in rows}})
        fr_meta.append({"id": "struct" if k == 0 else F["id"], "label": F["label"], "file": fname, "width": len(next(iter(F["gapped"].values()))), "split": len(next(iter(F["gapped"].values()))), "kind": F.get("kind", "sequence"), "short": F.get("short", F["id"])})
    meta = {"system": spec["system"], "pack": PACK, "homo": True, "title": spec["title"], "build": spec["build"], "n_rows": len(rows_out),
            "segments": {"A": {"name": spec["seg"], "sub": "A"}}, "frames": fr_meta,
            "orders": {"tree": order, "A": order, "default": "tree"}, "order_basis": spec.get("order_basis", "order of the alignment"), "order_short": spec.get("order_short", "tree"),
            "modules": {c: anchors for c in classes}, "module_ids": {str(v): k for k, v in MOD_ID.items()},
            "default": {"row": dflt["row"], "sub": "A", "col": dflt["col"], "why": dflt["why"]},
            "marks": spec.get("marks", []), "types": TYPES, "subs": SUBS, "networks": [k for k, v in sorted(NETS.items(), key=lambda kv: kv[1])],
            "kinds": ["within the toxin chain", "", "", "", "", ""], "contact_record": "<HHBBBBBH: i, j, kind|type<<3, sub, n, d10, net, chain mask>", "conf_codes": {"254": "experimental residue", "255": "not in structure"},
            "node_rule": spec["node_text"], "node_short": spec.get("node_short", "yes (contacts with another domain at a strict NODE column)"), "net_rule": "network = pair of modules (domains) joined by a contact",
            "colpairs": {"file": "colpairs.json", "min_support": minsup, "n_struct": dict(nstruct), "split": W}, "inputs": spec.get("inputs", {}), "cry": spec.get("cry", {}), "tangle": spec.get("tangle")}
    sizes["colpairs.json"] = dump("colpairs.json", {"split": W, "classes": {c: {"n_struct": nstruct[c], "min": minsup.get(c, 1), "pairs": colpairs.get(c, [])} for c in classes}})
    if spec.get("features"):
        sizes["features.json"] = dump("features.json", spec["features"])
    sizes["meta.json"] = dump("meta.json", meta); sizes["rows.json"] = dump("rows.json", rows_out)
    sizes["contacts.json"] = dump("contacts.json", contacts_out); sizes["residues.json"] = dump("residues.json", res_out)
    if verbose:
        print("written to", OUT)
        for k in sorted(sizes, key=lambda k: -sizes[k])[:8]:
            print("  %-22s %8.2f MB" % (k, sizes[k] / 1e6))
        print("  total %.2f MB in %d files; rows %d, with a structure %d, flagged %d, in %.0fs" % (sum(sizes.values()) / 1e6, len(sizes), len(rows_out), sum(1 for r in rows_out if r["bundle"]), sum(1 for r in rows_out if r["flags"]), time.time() - t0), flush=True)
    return {"res": res_out, "rows": rows_out, "cols": {r["id"]: r["cols"] for r in rows}, "strict": sorted(strict), "nets": NETS}
