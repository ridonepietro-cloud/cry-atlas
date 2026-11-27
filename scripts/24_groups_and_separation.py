"""Activity groups (nematode-active / insect-active ...) of the Cry rows and a test of which alignment columns and which residue-residue contacts separate nematode-active from insect-active toxins.
Groups (BPPRC specificity export, per row = union over its BPPRC names): 'Nematode-active' (an active record against Rhabditida), 'Insect-active' (active against an insect order, no nematode record),
'Other organism only', 'Tested, none active', 'No data'.
Signature test (structure-based frame columns with >= 30 % occupancy; contacts of the primary models):
  * a feature = (column, residue) or (column pair in contact);   positives = nematode-active rows of the six nematode families (Cry5, 12, 13, 14, 21, 95), negatives = insect-active rows;
  * candidate = present in >= 80 % of the positives, <= 10 % of the negatives, and in the majority of the rows of >= 5 of the 6 nematode families (the families are the units: relatedness would otherwise give 'markers' of one family);
  * validation = leave one NEMATODE FAMILY out: top-K features chosen on the other families and on half of the insect FAMILIES, scored (mean presence of the K features) on the held-out family against the other half of the insect families, AUC;
  * null = the same procedure with pseudo-nematode families drawn at random from the insect families (same number of families, similar number of rows): 'p' = share of random draws with a mean AUC >= the observed one.
Output: work/separation.json, docs/nematode_vs_insect_separation.md, docs/nematode_signature_columns.tsv.      Run: py 24_groups_and_separation.py     (numpy / scipy: not -I)"""
import os, sys, csv, json, collections, random, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import *
import cry_core as CC
t0 = time.time()
rows = jload(os.path.join(SITE, "data", "cry", "rows.json")); meta = jload(os.path.join(SITE, "data", "cry", "meta.json"))
fam = {r["id"]: (r.get("extra") or {}).get("family") for r in rows}; names_of = {r["id"]: (r.get("extra") or {}).get("names", []) for r in rows}
name2row = {n: i for i, ns in names_of.items() for n in ns}
INSECT = {"Lepidoptera", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Orthoptera", "Blattodea", "Siphonaptera", "Thysanoptera", "Isoptera", "Neuroptera", "Psocodea", "Odonata", "Mantodea", "Phasmatodea", "Dermaptera", "Ephemeroptera", "Trichoptera"}
act = collections.defaultdict(lambda: {"yes": collections.defaultdict(set), "no": collections.defaultdict(set)})
for x in csv.DictReader(open(os.path.join(CRY, "raw", "specificity", "association.csv"), encoding="utf-8", newline="")):
    i = name2row.get(x["name"]); a = (x["activity"] or "").strip().lower(); o = x["target_order"].strip()
    if i is None or a not in ("yes", "no") or not o: continue
    act[i][a][o].add(x["target_species"].strip())
group = {}; evid = {}
for r in rows:
    i = r["id"]; d = act.get(i)
    if not d: group[i] = "No data"; continue
    yes = set(d["yes"]);
    if "Rhabditida" in yes: group[i] = "Nematode-active"
    elif yes & INSECT: group[i] = "Insect-active"
    elif yes: group[i] = "Other organism only"
    else: group[i] = "Tested, none active"
    evid[i] = {"active_orders": {o: sorted(s) for o, s in d["yes"].items()}, "inactive_orders": {o: sorted(s) for o, s in d["no"].items()}}
cnt = collections.Counter(group.values()); print("activity groups:", dict(cnt))
NEMF = ["Cry5", "Cry12", "Cry13", "Cry14", "Cry21", "Cry95"]
pos_ids = [i for i in group if group[i] == "Nematode-active" and fam[i] in NEMF]; neg_ids = [i for i in group if group[i] == "Insect-active"]
print("positives (nematode families):", len(pos_ids), dict(collections.Counter(fam[i] for i in pos_ids)), "| negatives (insect-active):", len(neg_ids), "in", len({fam[i] for i in neg_ids}), "families")
# ---------- residue matrix on the structure frame
F0 = meta["frames"][0]; seqs = jload(os.path.join(SITE, "data", "cry", F0["file"]))["seqs"]; ids = [r["id"] for r in rows]; ix = {i: k for k, i in enumerate(ids)}
AA = "ACDEFGHIKLMNPQRSTVWY"; amap = np.full(256, 20, dtype=np.uint8)
for k, c in enumerate(AA): amap[ord(c)] = k
M = np.array([np.frombuffer(seqs[i].encode(), dtype=np.uint8) for i in ids]); M = amap[M]
occ = (M != 20).mean(0); cols = np.where(occ >= 0.30)[0]; Mc = M[:, cols]; print("frame columns", M.shape[1], "with occupancy >= 30 %:", len(cols))
from scipy import sparse
def onehot(Mc):
    n, c = Mc.shape; r = np.repeat(np.arange(n), c); cc = np.tile(np.arange(c), n); v = Mc.ravel()
    keep = v != 20
    return sparse.csr_matrix((np.ones(keep.sum(), dtype=np.float32), (r[keep], (cc * 20 + v)[keep])), shape=(n, c * 20))
XA = onehot(Mc)                                                          # rows x (columns * 20 residues)
# ---------- contacts (primary model of each row): pairs of frame columns
import base64, struct
con = jload(os.path.join(SITE, "data", "cry", "contacts.json")); resj = jload(os.path.join(SITE, "data", "cry", "residues.json"))
rowcols = {i: np.array(CC.col_index(seqs[i])) for i in ids}
pairidx = {}; rr = []; cc_ = []
has_model = []
for k, i in enumerate(ids):
    u = base64.b64decode(con[i]); n = len(u) // 11
    if n == 0: continue
    has_model.append(k); cs = rowcols[i]; seen = set()
    for q in range(n):
        a, b = struct.unpack_from("<HH", u, q * 11); p = (int(cs[a]), int(cs[b]))
        if p in seen: continue
        seen.add(p); j = pairidx.setdefault(p, len(pairidx)); rr.append(k); cc_.append(j)
XC = sparse.csr_matrix((np.ones(len(rr), dtype=np.float32), (rr, cc_)), shape=(len(ids), len(pairidx)))
freq_all = np.asarray(XC.sum(0)).ravel() / max(1, len(has_model)); keepc = np.where(freq_all >= 0.03)[0]; XC = XC[:, keepc]
pair_of = {j: p for p, j in pairidx.items()}; pairs_kept = [pair_of[j] for j in keepc]
print("rows with contacts", len(has_model), "| distinct column pairs", len(pairidx), "| with >= 3 % of the rows:", len(pairs_kept))
hm = set(has_model)

def select(X, pos, neg, fams_of_pos, top=None, strict=True):
    """feature indices ranked by (fraction of positives - fraction of negatives); strict: >= 80 % / <= 10 % and >= 5 of 6 (or all-but-one) positive families carry it in the majority of their rows"""
    fp = np.asarray(X[pos].mean(0)).ravel(); fn = np.asarray(X[neg].mean(0)).ravel(); sc = fp - fn
    ok = (fp >= 0.8) & (fn <= 0.1) if strict else np.ones_like(fp, bool)
    if strict and len(fams_of_pos) >= 2:
        need = max(1, len(fams_of_pos) - 1); cntf = np.zeros(X.shape[1], dtype=np.int32)
        for f, rws in fams_of_pos.items():
            cntf += (np.asarray(X[rws].mean(0)).ravel() >= 0.5).astype(np.int32)
        ok &= cntf >= need
    idx = np.where(ok)[0]; idx = idx[np.argsort(-sc[idx])]
    return idx if top is None else idx[:top], fp, fn

def auc(pos_s, neg_s):
    if len(pos_s) == 0 or len(neg_s) == 0: return float("nan")
    allv = np.concatenate([pos_s, neg_s]); rk = allv.argsort().argsort() + 1.0
    # average ranks for ties
    order = np.argsort(allv, kind="mergesort"); sv = allv[order]; ranks = np.empty(len(allv)); i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]: j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1; i = j + 1
    r = ranks[:len(pos_s)].sum(); return (r - len(pos_s) * (len(pos_s) + 1) / 2.0) / (len(pos_s) * len(neg_s))

rng = random.Random(20261127); K = 20
def lofo(X, pos_by_fam, neg_by_fam, usable=None, nsplit=4):
    """leave one positive family out; negatives split by family into train / test halves; mean AUC over folds (a fold needs >= 1 other positive family and >= 1 selected feature, else AUC 0.5)"""
    res = {}
    famlist = [f for f in pos_by_fam if len(pos_by_fam[f]) > 0]
    negfams = sorted(neg_by_fam)
    for f in famlist:
        train_pos = [i for g in famlist if g != f for i in pos_by_fam[g]]; fams_tr = {g: pos_by_fam[g] for g in famlist if g != f}; a_list = []
        for s in range(nsplit):
            rr_ = random.Random(1000 * s + len(f) + len(train_pos)); nf = negfams[:]; rr_.shuffle(nf); half = len(nf) // 2
            ntr = [i for g in nf[:half] for i in neg_by_fam[g]]; nte = [i for g in nf[half:] for i in neg_by_fam[g]]
            if not ntr or not nte: continue
            idx, _, _ = select(X, train_pos, ntr, fams_tr, top=K)
            if len(idx) == 0: a_list.append(0.5); continue
            Xs = X[:, idx]; sp = np.asarray(Xs[pos_by_fam[f]].mean(1)).ravel(); sn = np.asarray(Xs[nte].mean(1)).ravel(); a_list.append(auc(sp, sn))
        res[f] = float(np.nanmean(a_list)) if a_list else float("nan")
    return res

def run(X, valid_rows, name, nperm=150):
    """X: rows x features; valid_rows: row ids usable (a model for contacts)"""
    vr = set(valid_rows); P = [i for i in pos_ids if i in vr]; N = [i for i in neg_ids if i in vr]
    pk = [ix[i] for i in P]; nk = [ix[i] for i in N]
    pos_by_fam = collections.defaultdict(list); neg_by_fam = collections.defaultdict(list)
    for i in P: pos_by_fam[fam[i]].append(ix[i])
    for i in N: neg_by_fam[fam[i]].append(ix[i])
    idx_all, fp, fn = select(X, pk, nk, dict(pos_by_fam))
    obs = lofo(X, dict(pos_by_fam), dict(neg_by_fam)); mean_obs = float(np.nanmean(list(obs.values())))
    print("%s: positives %d (%d families), negatives %d (%d families) | candidates passing all criteria: %d | LOFO AUC per held-out family %s -> mean %.2f" % (name, len(P), len(pos_by_fam), len(N), len(neg_by_fam), len(idx_all), {k: round(v, 2) for k, v in obs.items()}, mean_obs), flush=True)
    # null: pseudo-nematode families drawn from the insect families (same number of families, similar number of rows)
    nfam = list(neg_by_fam); sizes = {f: len(v) for f, v in neg_by_fam.items()}; target_rows = len(P); nulls = []
    for p in range(nperm):
        pick = rng.sample(nfam, min(len(pos_by_fam), len(nfam))); tries = 0
        while sum(sizes[f] for f in pick) > 2.5 * target_rows and tries < 50: pick = rng.sample(nfam, min(len(pos_by_fam), len(nfam))); tries += 1
        pbf = {f: neg_by_fam[f] for f in pick}; nbf = {f: v for f, v in neg_by_fam.items() if f not in pick}
        if len(nbf) < 2: continue
        r_ = lofo(X, pbf, nbf, nsplit=2); nulls.append(float(np.nanmean(list(r_.values()))))
    nulls = np.array(nulls); pval = float((np.sum(nulls >= mean_obs) + 1) / (len(nulls) + 1))
    print("   null (%d pseudo-nematode draws): mean AUC %.2f, 95th percentile %.2f -> p = %.3f" % (len(nulls), nulls.mean(), np.percentile(nulls, 95), pval), flush=True)
    return {"positives": len(P), "negatives": len(N), "pos_families": {f: len(v) for f, v in pos_by_fam.items()}, "candidates": int(len(idx_all)), "auc_by_family": obs, "auc_mean": mean_obs, "null_mean": float(nulls.mean()), "null_p95": float(np.percentile(nulls, 95)), "p": pval, "n_null": int(len(nulls))}, idx_all, fp, fn

res_cols, idx_cols, fp_c, fn_c = run(XA, ids, "alignment columns (residue at column)")
res_con, idx_con, fp_k, fn_k = run(XC, [ids[k] for k in has_model], "contacts (column pairs)")
# ---------- tables
mods = {m["id"]: m["name"] for m in next(iter(meta["modules"].values()))}; colmod = {}
for m in next(iter(meta["modules"].values())):
    for c in range(m["cols"][0], m["cols"][1] + 1): colmod[c] = m["name"]
ref = next(r["id"] for r in rows if r["label"].startswith("Cry5Ba1")); ref1 = next(r["id"] for r in rows if r["label"].startswith("Cry1Aa1"))
def resnum(rid, col):                                                       # residue number (1-based) of the row at a frame column, '' when the row has a gap there
    cs = rowcols[rid]; w = np.where(cs == col)[0]; return int(w[0]) + 1 if len(w) else ""
out_cols = []
for f in idx_cols[:60]:
    c = int(cols[f // 20]); a = AA[f % 20]
    pf = {fm: float(np.mean([M[ix[i], c] == (f % 20) for i in pos_ids if fam[i] == fm])) for fm in NEMF}
    out_cols.append({"col": c + 1, "aa": a, "module": colmod.get(c + 1, "?"), "freq_pos": float(fp_c[f]), "freq_neg": float(fn_c[f]), "by_family": pf, "cry5Ba1_pos": resnum(ref, c + 1), "cry1Aa1_pos": resnum(ref1, c + 1)})
print("top candidate residues (column, residue, module, positives, negatives):", [(o["col"], o["aa"], o["module"], round(o["freq_pos"], 2), round(o["freq_neg"], 2)) for o in out_cols[:12]])
out_pairs = []
for f in idx_con[:60]:
    p = pairs_kept[f]; out_pairs.append({"cols": [p[0], p[1]], "modules": [colmod.get(p[0], "?"), colmod.get(p[1], "?")], "freq_pos": float(fp_k[f]), "freq_neg": float(fn_k[f]), "cry5Ba1_pos": [resnum(ref, p[0]), resnum(ref, p[1])]})
print("top candidate contacts:", [(o["cols"], o["modules"], round(o["freq_pos"], 2), round(o["freq_neg"], 2)) for o in out_pairs[:8]])
json.dump({"groups": group, "group_counts": dict(cnt), "evidence": evid, "positives": pos_ids, "negatives": neg_ids, "nematode_families": NEMF, "columns": res_cols, "contacts": res_con, "top_columns": out_cols, "top_contacts": out_pairs, "K": K,
           "method": __doc__}, open(os.path.join(WORK, "separation.json"), "w"), indent=0)
with open(os.path.join(CRY, "docs", "nematode_signature_columns.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["frame_col", "residue", "module", "freq_nematode", "freq_insect", "Cry5Ba1_residue_no", "Cry1Aa1_residue_no"] + ["freq_" + x for x in NEMF])
    for o in out_cols: w.writerow([o["col"], o["aa"], o["module"], round(o["freq_pos"], 3), round(o["freq_neg"], 3), o["cry5Ba1_pos"], o["cry1Aa1_pos"]] + [round(o["by_family"][x], 2) for x in NEMF])
print("done in %.0f s" % (time.time() - t0))
