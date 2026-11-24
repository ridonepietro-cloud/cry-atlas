"""Gene-synteny data of the Atlas genome map (Clinker style): global alignment identity between the proteins of the genome windows of all Atlas entries.

Input : site/data/fit/gc.json (loci with the annotated genes of every entry; atlas/build_gc_pack.py) and the GenBank windows of the curation tool (phylogeny_rebuild/gc_windows/<locus src>), which hold the translations.
Method: like Clinker (Gilchrist and Chooi 2021): every protein against every protein, global alignment (Needleman-Wunsch, BLOSUM62, gap open 10 / extend 1, Biopython PairwiseAligner = the parasail setting of Clinker),
        identity = identical positions / alignment columns (gaps included). All pairs are not aligned: all pairs are scored without traceback (global score) and only the pairs with a score >= 0.10 x the self score of the shorter protein are aligned for the identity (calibration on 3000 pairs: every pair with >= 30 % identity has a normalised score >= 0.16); identical sequences are 100 % without alignment; pairs with identity < 30 % are not stored (the Clinker default threshold).
Output: site/data/fit/gcsyn.json = { built, method, min_identity, ids: [gene key, ...], pairs: [i, j, identity, i, j, identity, ...] } with i < j indices into ids; gene key = protein id, or "<record>:<locus tag>" when the gene has none.
        analysis/genome_map/synteny_report.txt (counts). Run with py (about 10 minutes on 8 cores)."""
import os, sys, json, time, re, hashlib, collections
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PR = HERE
PACK = sys.argv[1] if len(sys.argv) > 1 else 'cry'
GCJ = os.path.join(HERE, 'site', 'data', PACK, 'gc.json'); OUT = os.path.join(HERE, 'site', 'data', PACK, 'gcsyn.json'); REP = os.path.join(HERE, 'docs', 'genome_map', 'synteny_report_%s.txt' % PACK)
os.makedirs(os.path.dirname(REP), exist_ok=True)
WIN = os.path.join(HERE, 'work', 'gc_windows')
MIN_ID = 30; NORM = 0.10
MURPHY = {}
for i, grp in enumerate(['LVIM', 'C', 'A', 'G', 'ST', 'P', 'FYW', 'EDNQ', 'KR', 'H']):
    for a in grp: MURPHY[a] = i

def gene_key(L, g):
    return g[3] if g[3] else ((L['acc'] + ':' + g[4]) if g[4] else None)

_STD = set('ACDEFGHIKLMNPQRSTVWYBZX')
def sanitize(s):                                                # rare letters of GenBank translations (U selenocysteine, O pyrrolysine, J, lower case) -> the BLOSUM62 alphabet
    return ''.join(c if c in _STD else {'U': 'C', 'O': 'K', 'J': 'L'}.get(c, 'X') for c in s.upper())

def load_sequences():
    from Bio import SeqIO
    gc = json.load(open(GCJ, encoding='utf-8')); seqs = {}; miss = 0; total = 0; cache = {}
    for eid, e in gc['entries'].items():
        for L in e.get('loci', []) or []:
            if L['src'] not in cache:
                byp = {}; byt = {}; bypos = {}; got = False
                for part in L['src'].split(' + '):                       # a window of two records (RefSeq and GenBank twins): the translations of both
                    p = os.path.join(WIN, part)
                    if not os.path.exists(p): continue
                    got = True; rec = SeqIO.read(p, 'genbank')
                    for f in rec.features:
                        if f.type != 'CDS' or 'translation' not in f.qualifiers: continue
                        tr = f.qualifiers['translation'][0].replace('*', '').strip()
                        for pid in f.qualifiers.get('protein_id', []): byp.setdefault(pid, tr)
                        for lt in f.qualifiers.get('locus_tag', []): byt.setdefault(lt, tr)
                        bypos.setdefault((int(f.location.start) + 1, int(f.location.end)), tr)
                cache[L['src']] = (byp, byt, bypos) if got else None
            c = cache[L['src']]
            for g in L['genes']:
                if g[8] & 4: continue                                  # RNA genes
                total += 1; k = gene_key(L, g)
                if k is None or k in seqs: continue
                s = None
                if c:
                    byp, byt, bypos = c
                    s = byp.get(g[3]) if g[3] else None
                    if s is None and g[4]: s = byt.get(g[4])
                    if s is None: s = bypos.get((g[0] - L['lo'] + 1, g[1] - L['lo'] + 1))
                if s: seqs[k] = sanitize(s)
                else: miss += 1
    return seqs, total, miss

_AL = None; _SEQ = None; _SELF = None
def _init(seqs=None):
    global _AL, _SEQ, _SELF
    from Bio.Align import PairwiseAligner, substitution_matrices
    _AL = PairwiseAligner(); _AL.mode = 'global'; _AL.substitution_matrix = substitution_matrices.load('BLOSUM62'); _AL.open_gap_score = -10; _AL.extend_gap_score = -1
    _SEQ = seqs
    if seqs is not None: _SELF = [_AL.score(x, x) for x in seqs]
def _scan(rng):                                                    # score-only pass for the rows i0 .. i1: pairs (i, j > i) whose global score is >= NORM x the self score of the shorter sequence
    i0, i1 = rng; out = []
    for i in range(i0, i1):
        si = _SEQ[i]; li = len(si)
        for j in range(i + 1, len(_SEQ)):
            sj = _SEQ[j]; lj = len(sj)
            if max(li, lj) > 3.2 * min(li, lj): continue            # identity over the alignment columns cannot reach 30 % with such a length ratio
            sc = _AL.score(si, sj)
            if sc >= NORM * min(_SELF[i], _SELF[j]): out.append((i, j))
    return out
def _ident(job):
    i, j, a, b = job
    al = _AL.align(a, b)[0]; (ca, cb) = al.aligned
    m = 0
    for (a0, a1), (b0, b1) in zip(ca, cb):
        m += sum(1 for x, y in zip(a[a0:a1], b[b0:b1]) if x == y)
    n = al.shape[1]
    return i, j, round(100.0 * m / n)

def main():
    t0 = time.time(); seqs, total, miss = load_sequences()
    keys = sorted(seqs); print('genes', total, 'with a sequence', len(keys), 'without', miss, file=sys.stderr)
    by_hash = collections.defaultdict(list)
    for i, k in enumerate(keys): by_hash[hashlib.md5(seqs[k].encode()).hexdigest()].append(i)
    reps = [v[0] for v in by_hash.values()]; print('unique sequences', len(reps), file=sys.stderr)
    seqL = [seqs[keys[r]] for r in reps]
    import multiprocessing as mp
    nproc = max(1, (os.cpu_count() or 4) - 1); n = len(seqL)
    # chunks of rows with about the same number of pairs
    cuts = [0]; tot = n * (n - 1) // 2; acc = 0; step = tot // (nproc * 12)
    for i in range(n):
        acc += n - i - 1
        if acc >= step: cuts.append(i + 1); acc = 0
    if cuts[-1] != n: cuts.append(n)
    rngs = list(zip(cuts[:-1], cuts[1:])); cand = []
    with mp.Pool(nproc, initializer=_init, initargs=(seqL,)) as pool:
        for k, r in enumerate(pool.imap_unordered(_scan, rngs)):
            cand += r; print('scan chunk', k + 1, 'of', len(rngs), 'candidates', len(cand), 'time', round(time.time() - t0), file=sys.stderr)
        jobs = [(a, b, seqL[a], seqL[b]) for a, b in cand]; res = []
        for k, r in enumerate(pool.imap_unordered(_ident, jobs, chunksize=100)):
            res.append(r)
            if k % 20000 == 0: print('identity', k, 'of', len(jobs), round(time.time() - t0), file=sys.stderr)
    members = {r: v for r, v in zip(reps, by_hash.values())}
    pairs = []; seen = set()
    def add(i, j, idn):
        if i == j: return
        a, b = (i, j) if i < j else (j, i)
        if (a, b) in seen: return
        seen.add((a, b)); pairs.append((a, b, idn))
    for v in by_hash.values():                                       # identical sequences
        for x in range(len(v)):
            for y in range(x + 1, len(v)): add(v[x], v[y], 100)
    ok = 0
    for a, b, idn in res:
        if idn < MIN_ID: continue
        ok += 1
        for i in members[reps[a]]:
            for j in members[reps[b]]: add(i, j, idn)
    pairs.sort(); flat = []
    for a, b, idn in pairs: flat += [a, b, idn]
    out = {'built': time.strftime('%Y-%m-%d'), 'method': 'global alignment (Needleman-Wunsch, BLOSUM62, gap open 10, extend 1), identity = identical positions / alignment columns; pairs are first scored (score only) and aligned for identity when the global score is >= %.2f x the self score of the shorter protein' % NORM, 'min_identity': MIN_ID, 'ids': keys, 'pairs': flat}
    json.dump(out, open(OUT, 'w'), separators=(',', ':'))
    hist = collections.Counter(min(9, p[2] // 10) for p in pairs)
    rep = ['genes (CDS) %d, with a sequence %d (%d without translation), unique sequences %d' % (total, len(keys), miss, len(reps)), 'candidate pairs aligned %d, pairs >= %d %% identity (incl. identical sequences, all gene keys) %d' % (len(jobs), MIN_ID, len(pairs)),
           'identity histogram (pairs): ' + ', '.join('%d-%d: %d' % (10 * b, 10 * b + 9, hist[b]) for b in sorted(hist)), 'file %s %d bytes, time %d s' % (OUT, os.path.getsize(OUT), time.time() - t0)]
    open(REP, 'w', encoding='utf-8').write('\n'.join(rep) + '\n'); print('\n'.join(rep))

if __name__ == '__main__':
    main()
