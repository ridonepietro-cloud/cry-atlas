"""Core helpers of the Cry Atlas pack builder (self-contained subset of atlas_common / homo_common / git_struct / fit_contacts / atlas_codec / build_fit of the NODE-MSA Atlas).
A ROW is the source of truth (ungapped sequence + per-residue data); an alignment frame is only a set of gapped strings; structure residues are placed on the row sequence (0-based index)."""
import os, re, json, base64, struct, collections
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "C", "PYL": "K"}
TYPES = ["HBOND", "VDW", "IONIC", "PIPISTACK", "PICATION", "PIHBOND", "SSBOND", "OTHER"]
SUBS = ["MC_MC", "MC_SC", "SC_MC", "SC_SC", "OTHER"]
FLIP = {"MC_SC": "SC_MC", "SC_MC": "MC_SC"}


def col_index(gapped):
    return [c + 1 for c, ch in enumerate(gapped) if ch != "-"]


def chain_sequence(residues):
    return "".join(AA3.get(r["resn"], "X") for r in residues)


def parse_pdb_chains(path):
    """ATOM records of a PDB file -> {chain: [ {resnum, icode, resn, atoms, bfac} ]} (first model, first altloc, heavy atoms)"""
    chains = {}
    for ln in open(path, encoding="utf-8", errors="replace"):
        if ln.startswith("ENDMDL"):
            break
        if not ln.startswith("ATOM"):
            continue
        if ln[16] not in " A":
            continue
        name = ln[12:16].strip()
        if name.startswith("H") or ln[76:78].strip() in ("H", "D"):
            continue
        ch = ln[21]; resnum = int(ln[22:26]); icode = ln[26].strip(); resn = ln[17:20]
        lst = chains.setdefault(ch, [])
        if not lst or lst[-1]["resnum"] != resnum or lst[-1]["icode"] != icode or lst[-1]["resn"] != resn:
            lst.append({"resnum": resnum, "icode": icode, "resn": resn, "atoms": {}, "bfac": float(ln[60:66])})
        lst[-1]["atoms"][name] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        if name == "CA":
            lst[-1]["bfac"] = float(ln[60:66])
    return chains


def parse_cif_chains(path, want=None):
    """mmCIF atom_site loop (heavy atoms, first model, first altloc); auth_asym_id / auth_seq_id, never label_*"""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    cols, i, n = [], 0, len(lines)
    while i < n:
        if lines[i].startswith("loop_") and i + 1 < n and lines[i + 1].startswith("_atom_site."):
            j = i + 1
            while j < n and lines[j].startswith("_atom_site."):
                cols.append(lines[j].strip()); j += 1
            i = j
            break
        i += 1
    else:
        raise ValueError("no atom_site loop in " + path)
    ix = {c: k for k, c in enumerate(cols)}
    for c in ["_atom_site.group_PDB", "_atom_site.label_atom_id", "_atom_site.label_comp_id", "_atom_site.auth_asym_id", "_atom_site.auth_seq_id", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
              "_atom_site.Cartn_z", "_atom_site.B_iso_or_equiv"]:
        if c not in ix:
            raise ValueError("missing column %s in %s" % (c, path))
    ic = ix.get("_atom_site.pdbx_PDB_ins_code"); al = ix.get("_atom_site.label_alt_id"); mo = ix.get("_atom_site.pdbx_PDB_model_num"); el = ix.get("_atom_site.type_symbol")
    chains, first_model = {}, None
    while i < n:
        ln = lines[i]; i += 1
        if ln.startswith("#") or ln.startswith("loop_") or ln.startswith("_"):
            break
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        t = ln.split()
        if len(t) < len(cols):
            continue
        if mo is not None:
            if first_model is None:
                first_model = t[mo]
            elif t[mo] != first_model:
                break
        if t[ix["_atom_site.group_PDB"]] != "ATOM":
            continue
        ch = t[ix["_atom_site.auth_asym_id"]]
        if want is not None and ch not in want:
            continue
        if al is not None and t[al] not in (".", "?", "A"):
            continue
        if (t[el] if el is not None else "") in ("H", "D"):
            continue
        name = t[ix["_atom_site.label_atom_id"]].strip('"'); resn = t[ix["_atom_site.label_comp_id"]]; resnum = int(t[ix["_atom_site.auth_seq_id"]])
        icode = "" if ic is None or t[ic] in (".", "?") else t[ic]
        lst = chains.setdefault(ch, [])
        if not lst or lst[-1]["resnum"] != resnum or lst[-1]["icode"] != icode or lst[-1]["resn"] != resn:
            lst.append({"resnum": resnum, "icode": icode, "resn": resn, "atoms": {}, "bfac": float(t[ix["_atom_site.B_iso_or_equiv"]])})
        lst[-1]["atoms"][name] = (float(t[ix["_atom_site.Cartn_x"]]), float(t[ix["_atom_site.Cartn_y"]]), float(t[ix["_atom_site.Cartn_z"]]))
        if name == "CA":
            lst[-1]["bfac"] = float(t[ix["_atom_site.B_iso_or_equiv"]])
    return chains


def read_chains(path, want=None):
    ch = parse_pdb_chains(path) if path.lower().endswith(".pdb") else parse_cif_chains(path, want)
    return {c: v for c, v in ch.items() if want is None or c in want}


_aligner = None


def map_chain_to_seq(chain_seq, seq):
    """residues of a structure chain (file order) -> 0-based indices of seq: (idx, identity over placed, mode)"""
    global _aligner
    if not chain_seq or not seq:
        return [-1] * len(chain_seq), 0.0, "empty"
    at = seq.find(chain_seq)
    if at >= 0:
        return list(range(at, at + len(chain_seq))), 1.0, "identical stretch"
    if _aligner is None:
        from Bio.Align import PairwiseAligner
        _aligner = PairwiseAligner(); _aligner.mode = "global"; _aligner.match_score, _aligner.mismatch_score = 2, -1
        _aligner.open_gap_score, _aligner.extend_gap_score = -6, -1; _aligner.end_gap_score = 0
    aln = _aligner.align(seq, chain_seq)[0]
    idx, same, tot = [-1] * len(chain_seq), 0, 0
    for (s0, s1), (c0, c1) in zip(*aln.aligned):
        s0, s1, c0, c1 = int(s0), int(s1), int(c0), int(c1)
        for k in range(s1 - s0):
            idx[c0 + k] = s0 + k; tot += 1; same += seq[s0 + k] == chain_seq[c0 + k]
    return idx, (same / tot if tot else 0.0), "alignment"


def cif_scheme(path):
    """entities, SEQRES (canonical one-letter), UniProt refs and the residue scheme (auth_seq_id + icode -> 0-based SEQRES position) of an mmCIF file"""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    d = MMCIF2Dict(path)
    g = lambda k: (lambda v: v if isinstance(v, list) else [v])(d.get(k, []))
    ent_seq = {e: s.replace("\n", "").replace(" ", "") for e, s in zip(g("_entity_poly.entity_id"), g("_entity_poly.pdbx_seq_one_letter_code_can"))}
    ent_name = dict(zip(g("_entity.id"), g("_entity.pdbx_description")))
    ref = {}
    for e, db, acc in zip(g("_struct_ref.entity_id"), g("_struct_ref.db_name"), g("_struct_ref.pdbx_db_accession")):
        ref.setdefault(e, []).append((db, acc))
    chains = {}
    for asym, ent, seqid, authnum, ins, authchain in zip(g("_pdbx_poly_seq_scheme.asym_id"), g("_pdbx_poly_seq_scheme.entity_id"), g("_pdbx_poly_seq_scheme.seq_id"),
                                                         g("_pdbx_poly_seq_scheme.pdb_seq_num"), g("_pdbx_poly_seq_scheme.pdb_ins_code"), g("_pdbx_poly_seq_scheme.pdb_strand_id")):
        c = chains.setdefault(authchain, {"entity": ent, "asym": asym, "map": {}})
        if authnum in ("?", "."):
            continue
        try:
            c["map"][(int(authnum), "" if ins in ("?", ".") else ins)] = int(seqid) - 1
        except ValueError:
            pass
    title = (g("_struct.title") or [""])[0]
    return {"entities": {e: {"seq": ent_seq.get(e, ""), "name": ent_name.get(e, ""), "ref": ref.get(e, [])} for e in ent_seq}, "chains": chains, "title": title}


def _ident(res, idx, seq):
    cs = chain_sequence(res); same = sum(1 for k, i in enumerate(idx) if i >= 0 and seq[i] == cs[k]); tot = sum(1 for i in idx if i >= 0)
    return (same / tot if tot else 0.0), tot


def place_chain(res, seq, chain=None, scheme=None):
    """0-based indices of the residues of a chain on `seq`: (1) exactly through the SEQRES scheme when the entity is `seq`; (2) author numbering when residues agree; (3) global alignment"""
    if scheme and chain in scheme["chains"]:
        e = scheme["entities"].get(scheme["chains"][chain]["entity"])
        if e and e["seq"] == seq:
            m = scheme["chains"][chain]["map"]; idx = [m.get((x["resnum"], x["icode"]), -1) for x in res]
            ident, tot = _ident(res, idx, seq)
            if tot and ident >= 0.97:
                return idx, ident, "SEQRES scheme"
    idx = [x["resnum"] - 1 if not x["icode"] and 0 < x["resnum"] <= len(seq) else -1 for x in res]
    ident, tot = _ident(res, idx, seq)
    if tot >= 0.9 * len(res) and ident >= 0.97:
        return idx, ident, "author numbering"
    return map_chain_to_seq(chain_sequence(res), seq)


def one_chain_structure(path, chain, seq, scheme=None):
    """the protomer `chain` of a structure placed on the row sequence (chain dict of the Atlas)"""
    ch = read_chains(path, {chain}); res = ch[chain]; idx, ident, mode = place_chain(res, seq, chain, scheme)
    c = {"sub": "A", "res": res, "seq": chain_sequence(res), "idx": idx, "ident": ident, "mode": mode, "placed": sum(1 for i in idx if i >= 0),
         "map": {(x["resnum"], x["icode"]): i for x, i in zip(res, idx) if i >= 0}}
    return {"chains": {chain: c}, "sub_of": {chain: "A"}, "order": {"A": [chain], "B": []}, "rep_A": chain, "rep_B": None}


def edges_intra(path, S):
    """RING edge file -> {(i, j, main, sub): [dmin, count]} for contacts INSIDE the protomer chain of S (i < j; MC_SC / SC_MC flipped to keep i < j); stats"""
    ch = S["order"]["A"][0]; mp = S["chains"][ch]["map"]; acc, st = {}, collections.Counter()
    with open(path, encoding="utf-8") as fh:
        next(fh)
        for ln in fh:
            t = ln.rstrip("\n").split("\t")
            if len(t) < 4:
                continue
            st["edges"] += 1
            c1, r1, ic1, _ = t[0].split(":"); c2, r2, ic2, _ = t[2].split(":")
            if c1 != ch or c2 != ch:
                st["other_chain"] += 1; continue
            i1 = mp.get((int(r1), "" if ic1 == "_" else ic1)); i2 = mp.get((int(r2), "" if ic2 == "_" else ic2))
            if i1 is None or i2 is None:
                st["unplaced_residue"] += 1; continue
            main, _, sub = t[1].partition(":")
            if main not in TYPES:
                main = "OTHER"
            if sub not in SUBS:
                sub = "OTHER"
            if i1 > i2:
                i1, i2, sub = i2, i1, FLIP.get(sub, sub)
            if i1 == i2:
                continue
            key = (i1, i2, main, sub); a = acc.setdefault(key, [float(t[3]), 0]); a[0] = min(a[0], float(t[3])); a[1] += 1; st["kept"] += 1
    return acc, st


# ---- pack encodings (as build_fit.py / atlas_codec.py) ----
ATOMS = {"A": "N CA C O CB", "R": "N CA C O CB CG CD NE CZ NH1 NH2", "N": "N CA C O CB CG OD1 ND2", "D": "N CA C O CB CG OD1 OD2", "C": "N CA C O CB SG", "Q": "N CA C O CB CG CD OE1 NE2",
         "E": "N CA C O CB CG CD OE1 OE2", "G": "N CA C O", "H": "N CA C O CB CG ND1 CD2 CE1 NE2", "I": "N CA C O CB CG1 CG2 CD1", "L": "N CA C O CB CG CD1 CD2", "K": "N CA C O CB CG CD CE NZ",
         "M": "N CA C O CB CG SD CE", "F": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ", "P": "N CA C O CB CG CD", "S": "N CA C O CB OG", "T": "N CA C O CB OG1 CG2",
         "W": "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2", "Y": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH", "V": "N CA C O CB CG1 CG2"}
ATOMS = {k: v.split() for k, v in ATOMS.items()}; GENERIC = ["N", "CA", "C", "O", "CB"]; ALIAS = {"SE": "SD"}


def atom_list(aa):
    return ATOMS.get(aa, GENERIC)


def enc_chain(seq, residues, bfac):
    out, prev, masks = bytearray(), None, {}
    for pos, (aa, atoms) in enumerate(zip(seq, residues)):
        names = atom_list(aa); have = {ALIAS.get(k, k): v for k, v in atoms.items()}; bits = 0
        for bit, nm in enumerate(names):
            at = have.get(nm)
            if at is None:
                continue
            bits |= 1 << bit; q = [int(round(v * 10)) for v in at]
            if prev is not None and all(-127 <= q[k] - prev[k] <= 127 for k in range(3)):
                out.extend(struct.pack("<3b", *(q[k] - prev[k] for k in range(3))))
            else:
                if any(abs(v) > 32767 for v in q):
                    raise ValueError("coordinate out of int16 range")
                out.append(0x80); out.extend(struct.pack("<3h", *q))
            prev = q
        if bits != (1 << len(names)) - 1:
            masks[pos] = bits
    rec = {"enc": base64.b64encode(bytes(out)).decode("ascii"), "b": base64.b64encode(bytes(max(0, min(255, int(round(x)))) for x in bfac)).decode("ascii")}
    if masks:
        rec["m"] = masks
    return rec


def dec_chain(seq, rec):
    u = base64.b64decode(rec["enc"]); pts, i, p = [], 0, [0, 0, 0]
    while i < len(u):
        if u[i] == 0x80:
            p = list(struct.unpack_from("<3h", u, i + 1)); i += 7
        else:
            d = struct.unpack_from("<3b", u, i); p = [p[k] + d[k] for k in range(3)]; i += 3
        pts.append((p[0] / 10, p[1] / 10, p[2] / 10))
    masks, rows, j = rec.get("m") or {}, [], 0
    for pos, aa in enumerate(seq):
        names = atom_list(aa); bits = masks.get(pos, masks.get(str(pos), (1 << len(names)) - 1)); d = {}
        for bit, nm in enumerate(names):
            if bits >> bit & 1:
                d[nm] = pts[j]; j += 1
        rows.append(d)
    assert j == len(pts), "atom count mismatch"
    return rows


def b64(b):
    return base64.b64encode(bytes(b)).decode("ascii")


def pack_contacts(recs):
    out = bytearray()
    for i, j, kind, typ, sub, n, d10, net, cm in sorted(recs, key=lambda r: (r[2], r[0], r[1], r[3])):
        out += struct.pack("<HHBBBBBH", i, j, kind | (typ << 3), sub, min(n, 255), min(d10, 255), net, cm)
    return b64(out)


def runs_of(pairs):
    out = []
    for p, v in pairs:
        if out and p == out[-1][0] + out[-1][2] and v == out[-1][1] + out[-1][2]:
            out[-1][2] += 1
        else:
            out.append([p, v, 1])
    return out


def chain_records(S, chs):
    out = []
    for ch in chs:
        c = S["chains"][ch]; aa = chain_sequence(c["res"])
        rec = enc_chain(aa, [x["atoms"] for x in c["res"]], [x["bfac"] for x in c["res"]])
        rec.update({"seq": aa, "chain": ch, "ident": round(c["ident"], 4), "mode": c["mode"], "placed": c["placed"],
                    "rn": runs_of(list(enumerate(x["resnum"] for x in c["res"]))), "ix": runs_of([(p, i) for p, i in enumerate(c["idx"]) if i >= 0])})
        out.append(rec)
    return out


def read_fasta(path):
    out, h, s = [], None, []
    for ln in open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if h is not None:
                out.append((h, "".join(s)))
            h, s = ln[1:], []
        elif ln:
            s.append(ln.strip())
    if h is not None:
        out.append((h, "".join(s)))
    return out
