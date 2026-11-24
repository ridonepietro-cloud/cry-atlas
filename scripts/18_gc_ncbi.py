"""Cry Atlas copy of build_gc_ncbi.py (NODE-MSA Atlas): the genomic context of every row straight from NCBI.   Run: py -I 18_gc_ncbi.py cry   (input work/cry/gc_targets.json)

Input : work/<pack>/gc_targets.json  { <row id>: { "accs": [protein accessions: RefSeq WP_/NP_/XP_, GenBank, UniProt], "org": "species of the entry", "label": "..." } }  (written by the pack builder)
Method: efetch protein rettype=ipg (identical protein group) lists the genome records that encode the protein (nucleotide accession, start, stop, strand, organism, assembly); the record is chosen by species
        match to the entry, then RefSeq NC_ over NZ_ over GenBank; the window = gene +- PAD (15 kb, as the curation tool) is fetched from nuccore (efetch rettype=gb, seq_start / seq_stop) into
        phylogeny_rebuild/gc_windows/<acc>_<lo>_<hi>.gb (the folder build_gc_syn.py reads the translations from), and parsed like build_genome_context.py of the curation tool.
Output: site/data/<pack>/gc.json = the format of build_gc_pack.py: { built, source, pack, pad_bp, entries: { <row id>: { item, a, b, type, how, loci: [{acc, org, lo, hi, rec_lo, rec_hi, genes, roles: {a: {...}}, src, fd}], why } } }
        analysis/genome_map/gc_ncbi_report_<pack>.tsv.   IPG answers are cached in work/<pack>/ipg/.   Run with py.   NCBI allows 3 requests per second without a key."""
import os, re, sys, json, time, datetime, threading, urllib.request, urllib.parse, urllib.error
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PR = HERE
PACK = sys.argv[1] if len(sys.argv) > 1 else 'cry'; PAD = 15000
WORK = os.path.join(HERE, 'work', PACK); IPGD = os.path.join(WORK, 'ipg'); os.makedirs(IPGD, exist_ok=True)
GCW = os.path.join(HERE, 'work', 'gc_windows'); os.makedirs(GCW, exist_ok=True)
OUT = os.environ.get('GC_OUT') or os.path.join(HERE, 'site', 'data', PACK, 'gc.json'); REP = os.path.join(HERE, 'docs', 'genome_map'); os.makedirs(REP, exist_ok=True); os.makedirs(os.path.dirname(OUT), exist_ok=True)
EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
lock = threading.Lock(); last = [0.0]
def get(url, tries=5, timeout=180):
    err = None
    for i in range(tries):
        with lock:
            w = 0.36 - (time.time() - last[0])
            if w > 0: time.sleep(w)
            last[0] = time.time()
        try: return urllib.request.urlopen(url, timeout=timeout).read().decode('utf-8', 'replace')
        except urllib.error.HTTPError as e:
            err = e
            if e.code in (400, 404): break
            time.sleep(4 * (i + 1))
        except Exception as e: err = e; time.sleep(4 * (i + 1))
    raise err
def ipg(acc):
    fn = os.path.join(IPGD, re.sub(r'[^A-Za-z0-9_.]', '_', acc) + '.tsv')
    if os.path.exists(fn) and os.path.getsize(fn) > 20: txt = open(fn, encoding='utf-8').read()
    else:
        try: txt = get(EU + urllib.parse.urlencode({'db': 'protein', 'id': acc, 'rettype': 'ipg', 'retmode': 'text'}))
        except Exception: txt = ''
        open(fn, 'w', encoding='utf-8').write(txt)
    rows = []
    for ln in txt.splitlines()[1:]:
        p = ln.split('\t')
        if len(p) < 9 or not p[3].isdigit(): continue
        rows.append(dict(src=p[1], nuc=p[2], s=int(p[3]), e=int(p[4]), strand=p[5], prot=p[6], name=p[7], org=p[8], strain=p[9] if len(p) > 9 else '', asm=p[10] if len(p) > 10 else ''))
    return rows
def words(s): return set(w for w in re.findall(r'[a-z0-9]+', (s or '').lower()) if len(w) > 2 and w not in ('str', 'strain', 'substr', 'subsp', 'bacterium', 'uncultured'))
def pick(rows, org):
    ow = words(org)
    def score(r):
        rw = words(r['org'] + ' ' + r['strain']); ov = len(ow & rw) / max(1, len(ow)) if ow else 0
        return (-round(ov, 2), 0 if r['nuc'].startswith('NC_') else 1 if r['nuc'].startswith('NZ_') else 2 if r['src'] == 'RefSeq' else 3)
    return sorted(rows, key=score)[0] if rows else None
# GenBank windows, parsed as the curation tool does
from Bio import SeqIO
def fetch_window(acc, lo, hi):
    fn = os.path.join(GCW, '%s_%d_%d.gb' % (acc, lo, hi))
    if os.path.exists(fn) and os.path.getsize(fn) > 500: return fn
    t = get(EU + urllib.parse.urlencode({'db': 'nuccore', 'id': acc, 'rettype': 'gb', 'retmode': 'text', 'seq_start': lo, 'seq_stop': hi}))
    if not t.startswith('LOCUS'): raise RuntimeError('no GenBank text: ' + t[:80])
    open(fn, 'w', encoding='utf-8').write(t); return fn
def parse_window(fn):
    rec = SeqIO.read(fn, 'genbank'); m = re.search(r'REGION: (\d+)\.\.(\d+)', open(fn, encoding='utf-8').read(2000)); lo = int(m.group(1)) if m else 1; off = lo - 1
    genes = []; trk = {}
    for f in rec.features:
        if f.type not in ('CDS', 'tRNA', 'rRNA', 'ncRNA', 'tmRNA'): continue
        q = f.qualifiers; s = int(f.location.start) + 1 + off; e = int(f.location.end) + off; st = 1 if (f.location.strand or 1) > 0 else -1
        tr = (q.get('translation') or [''])[0]; aa = len(tr) if tr else (e - s + 1) // 3 - 1
        pseudo = ('pseudo' in q) or ('pseudogene' in q); part = ('<' in str(f.location)) or ('>' in str(f.location))
        trk[(s, e)] = tr
        genes.append([s, e, st, (q.get('protein_id') or [''])[0], (q.get('locus_tag') or [''])[0], (q.get('gene') or [''])[0], ((q.get('product') or [''])[0] or (f.type if f.type != 'CDS' else ''))[:90], aa if f.type == 'CDS' else 0, (1 if pseudo else 0) | (2 if part else 0) | (4 if f.type != 'CDS' else 0)])
    genes.sort(key=lambda g: (g[0], g[1]))
    return dict(acc=rec.id, lo=lo, hi=lo + len(rec.seq) - 1, org=rec.annotations.get('organism', ''), genes=genes, trk=trk, file=os.path.basename(fn), fd=datetime.date.fromtimestamp(os.path.getmtime(fn)).isoformat())
def find_gene(genes, loc):
    best = None
    for i, g in enumerate(genes):
        ov = min(g[1], loc['e']) - max(g[0], loc['s']) + 1
        if ov <= 0: continue
        sc = ov / max(1, max(g[1] - g[0] + 1, loc['e'] - loc['s'] + 1)) + (0.2 if (g[2] > 0) == (loc['strand'] == '+') else 0) + (0.5 if loc.get('pid') and g[3] == loc['pid'] else 0) - (0.3 if g[8] & 4 else 0)
        if sc > 0.5 and (best is None or sc > best[0]): best = (sc, i)
    return best[1] if best else None
def pdbe_uniprot(pdb):
    fn = os.path.join(IPGD, 'pdb_' + pdb + '.json')
    if os.path.exists(fn): return json.load(open(fn))
    try: out = list(json.loads(get('https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/%s' % pdb.lower())).get(pdb.lower(), {}).get('UniProt', {}).keys())
    except Exception: out = []
    json.dump(out, open(fn, 'w')); return out
MISSING = False  # (curation-tool merge mode of the original script removed)                              # --missing: only the rows of the pack that have no genome map yet; the result is merged into the existing gc.json (curation-tool windows are kept)
if MISSING:
    rows_ = json.load(open(os.path.join(HERE, 'site', 'data', PACK, 'rows.json'), encoding='utf-8')); old_ = json.load(open(OUT, encoding='utf-8'))['entries'] if os.path.exists(OUT) else {}
    targets = {}
    _meta = json.load(open(os.path.join(HERE, 'site', 'data', PACK, 'meta.json'), encoding='utf-8')); _aln = json.load(open(os.path.join(HERE, 'site', 'data', PACK, 'aln_structural.json'), encoding='utf-8'))['seqs']
    _sp = _meta['frames'][0]['split']; _w = _meta['frames'][0]['width']; SEQS_ = {}
    for rid_, g_ in _aln.items():
        a_, b_ = g_[:_sp].replace('-', ''), g_[_sp:].replace('-', '')
        SEQS_[rid_] = {'seqA': a_, **({'seqB': b_} if b_ and _sp < _w else {})}
    for r_ in rows_:
        if (old_.get(r_['id']) or {}).get('loci') and (old_[r_['id']].get('how') != 'NCBI identical-protein group' or not SEQS_.get(r_['id'], {}).get('seqB')): continue       # rows done by NCBI are redone when the entry has a B subunit (to mark the B gene)
        accs_ = [a for a in (r_.get('acc'), r_.get('accB')) if a and re.match(r'^[A-Z]{1,3}_?\d{5,}(\.\d+)?$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$', a)]
        ov_ = (json.load(open(os.path.join(HERE, 'entries', 'gc_acc_overrides.json'), encoding='utf-8')).get(PACK) or {}) if os.path.exists(os.path.join(HERE, 'entries', 'gc_acc_overrides.json')) else {}
        if r_['id'] in ov_ and isinstance(ov_[r_['id']], list): accs_ = list(ov_[r_['id']])
        elif r_['id'] in ov_ and isinstance(ov_[r_['id']], dict): targets[r_['id']] = {'accs': [], 'org': r_.get('species') or '', 'label': r_['label'], 'locus': ov_[r_['id']]['locus']}; continue
        if not accs_ and r_.get('pdb'): accs_ = pdbe_uniprot(r_['pdb'])[:2]
        accs_ = list(dict.fromkeys(accs_)); targets[r_['id']] = {'accs': accs_, 'org': r_.get('species') or '', 'label': r_['label']}
        if SEQS_.get(r_['id']): targets[r_['id']].update(SEQS_[r_['id']])
else:
    targets = json.load(open(os.path.join(WORK, 'gc_targets.json'), encoding='utf-8'))
UNIP = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
def uniprot_refseq(acc):
    """RefSeq protein accessions of a UniProt entry (cached): the identical protein group of a UniProt accession can list only mRNA / cDNA records"""
    fn = os.path.join(IPGD, acc + '.uniprot2.json')
    if os.path.exists(fn): return json.load(open(fn))
    try:
        d = json.loads(get('https://rest.uniprot.org/uniprotkb/%s.json' % acc)); xr = d.get('uniProtKBCrossReferences', [])
        out = [x['id'] for x in xr if x['database'] == 'RefSeq']
        for x in xr:                                                      # INSDC protein ids of the EMBL cross-references (TrEMBL entries often have no RefSeq cross-reference)
            if x['database'] == 'EMBL':
                for p_ in x.get('properties', []):
                    if p_.get('key') == 'ProteinId' and p_.get('value') not in ('-', None) and p_['value'] not in out: out.append(p_['value'])
    except Exception: out = []
    json.dump(out, open(fn, 'w')); return out
def uniparc_ids(acc):
    """protein ids (RefSeq / INSDC) of the UniParc entry of a UniProt accession: the identical sequences of other records, often with a genome (cached)"""
    fn = os.path.join(IPGD, acc + '.uniparc.json')
    if os.path.exists(fn): return json.load(open(fn))
    out = []
    try:
        d = json.loads(get('https://rest.uniprot.org/uniprotkb/%s.json' % acc)); upi = (d.get('extraAttributes') or {}).get('uniParcId')
        if upi:
            u = json.loads(get('https://rest.uniprot.org/uniparc/%s.json' % upi))
            for x in u.get('uniParcCrossReferences', []):
                if x.get('database') in ('RefSeq', 'EMBL', 'EMBLWGS') and x.get('active', True) and x.get('id') and x['id'] not in out: out.append(x['id'])
    except Exception: pass
    out = out[:14]; json.dump(out, open(fn, 'w')); return out
def candidates(accs):
    """IPG records of the accessions (and, for UniProt accessions, of their RefSeq proteins), best first"""
    seen = []; allrows = []
    for a in accs:
        for q in [a] + (uniprot_refseq(a) if UNIP.match(a) else []):
            if q in seen: continue
            seen.append(q)
            for r in ipg(q): allrows.append((q, r))
    if not allrows:                                                  # UniProt accession without a genome record in its own identical protein group: the identical sequences of UniParc
        for a in accs:
            if UNIP.match(a):
                for q in uniparc_ids(a):
                    if q in seen: continue
                    seen.append(q)
                    for r in ipg(q): allrows.append((q, r))
                    if len(allrows) > 3: break
    return allrows
_AL = None
def ident(a, b):
    global _AL
    from Bio import Align
    from Bio.Align import substitution_matrices
    if _AL is None:
        _AL = Align.PairwiseAligner(); _AL.substitution_matrix = substitution_matrices.load('BLOSUM62'); _AL.open_gap_score = -10; _AL.extend_gap_score = -0.5; _AL.mode = 'global'
    if not a or not b or max(len(a), len(b)) > 3.0 * min(len(a), len(b)): return 0
    x = _AL.align(a, b)[0]; same = 0
    for (s0, s1), (t0, t1) in zip(*x.aligned):
        for k in range(s1 - s0): same += a[s0 + k] == b[t0 + k]
    return round(100 * same / max(len(a), len(b)))
def one(item):
    rid, t = item
    if t.get('locus'):                                                  # explicit location (accession, start, stop, strand, note): used when the protein has no genome record any more
        nuc, s0, e0, st0, note = t['locus']; lo = max(1, s0 - PAD); hi = e0 + PAD
        try: fn = fetch_window(nuc, lo, hi)
        except Exception as ex: return rid, {'why': 'The annotation of %s could not be retrieved (%s).' % (nuc, str(ex)[:60])}, 'fetch failed'
        w = parse_window(fn); lo2 = max(w['lo'], lo); hi2 = min(w['hi'], hi); genes = [g for g in w['genes'] if g[1] >= lo2 and g[0] <= hi2]
        x = dict(s=s0, e=e0, strand=st0, pid=None, how=note); i = find_gene(genes, x)
        L = dict(acc=w['acc'], org=w['org'], lo=lo2, hi=hi2, rec_lo=w['lo'], rec_hi=w['hi'], genes=genes, roles={'a': dict(i=i, s=s0, e=e0, strand=st0, pid=(genes[i][3] if i is not None else None), how=note)}, src=w['file'], fd=w['fd'])
        return rid, {'item': None, 'a': None, 'b': None, 'type': 'single', 'how': 'gene location given by hand', 'loci': [L], 'notes': [note]}, '%s:%d-%d (%d genes)' % (w['acc'], s0, e0, len(genes))
    cands = candidates(t.get('accs', []))
    if not cands: return rid, {'why': 'NCBI lists no genome record for the protein (%s).' % ', '.join(t.get('accs', [])) if t.get('accs') else 'The entry has no protein accession.'}, 'no IPG'
    ow = words(t.get('org') or '')
    def score(c):
        r = c[1]; rw = words(r['org'] + ' ' + r['strain']); ov = len(ow & rw) / max(1, len(ow)) if ow else 0
        return (-round(ov, 2), 0 if r['nuc'].startswith('NC_') else 1 if r['nuc'].startswith('NZ_') else 2 if r['src'] == 'RefSeq' else 3)
    cands.sort(key=score); best = None; tried = 0
    for used, r in cands:
        if tried >= 6: break
        tried += 1
        lo = max(1, min(r['s'], r['e']) - PAD); hi = max(r['s'], r['e']) + PAD
        try: fn = fetch_window(r['nuc'], lo, hi)
        except Exception: continue
        w = parse_window(fn); ctx = (w['hi'] - w['lo'] + 1) - (abs(r['e'] - r['s']) + 1)
        if best is None or ctx > best[0]: best = (ctx, used, r, w, lo, hi)
        if ctx >= 5000: break                                          # a genomic record with context on the sides (a short mRNA / cDNA record has none)
    if best is None: return rid, {'why': 'The annotation of the genome records of the protein could not be retrieved.'}, 'fetch failed'
    ctx, used, r, w, lo, hi = best
    lo2 = max(w['lo'], lo); hi2 = min(w['hi'], hi); genes = [g for g in w['genes'] if g[1] >= lo2 and g[0] <= hi2]
    x = dict(s=min(r['s'], r['e']), e=max(r['s'], r['e']), strand=r['strand'], pid=r['prot'], how='NCBI identical-protein group of ' + str(used))
    i = find_gene(genes, x)
    L = dict(acc=w['acc'], org=w['org'], lo=lo2, hi=hi2, rec_lo=w['lo'], rec_hi=w['hi'], genes=genes, roles={'a': dict(i=i, s=x['s'], e=x['e'], strand=x['strand'], pid=x['pid'], how=x['how'])}, src=w['file'], fd=w['fd'])
    bnote = None
    if t.get('seqB') and i is not None:                                   # the B gene: the CDS near the A gene whose translation is the B subunit of the row (sequence match; the accession of the B is not stored for these rows)
        best = None
        for j in range(max(0, i - 6), min(len(genes), i + 7)):
            g = genes[j]; tj = w['trk'].get((g[0], g[1]))
            if j == i or not tj or g[8] & 4 or g[2] != genes[i][2]: continue
            idn = ident(t['seqB'], tj.replace('*', ''))
            if idn >= 40 and (best is None or idn > best[0]): best = (idn, j)
        if best:
            g = genes[best[1]]; L['roles']['b'] = dict(i=best[1], s=g[0], e=g[1], strand='+' if g[2] > 0 else '-', pid=g[3] or g[4], how='sequence match to the B subunit of the entry (%d %% identical)' % best[0])
        else: bnote = 'No CDS next to the A gene matches the B subunit of the entry (identity >= 40 %): the B gene is not marked.'
    if t.get('seqA') and i is not None:
        ta = w['trk'].get((genes[i][0], genes[i][1])); ia = ident(t['seqA'], ta.replace('*', '')) if ta else 0
        if ia < 90: bnote = ((bnote + ' ') if bnote else '') + 'The A gene of the window is %d %% identical to the A subunit of the entry.' % ia
    e = {'item': None, 'a': used, 'b': (L['roles'].get('b') or {}).get('pid'), 'type': 'pair' if 'b' in L['roles'] else 'single', 'how': 'NCBI identical-protein group', 'loci': [L], 'notes': []}
    if bnote: e['notes'].append(bnote)
    if ctx < 5000: e['notes'].append('Only a short record (%s, %d bp) is listed for this protein: the window has no genomic context.' % (w['acc'], w['hi'] - w['lo'] + 1))
    if r['org'] and t.get('org') and not (words(t['org']) & words(r['org'])): e['notes'].append('Genome of %s (the entry is %s): the identical protein group has no record of the entry species.' % (r['org'], t['org']))
    return rid, e, '%s:%d-%d %s (%d genes)%s' % (w['acc'], x['s'], x['e'], r['org'], len(genes), '' if i is not None else ' NO CDS AT THE LOCUS')
entries = {}; rep = ['row\tlabel\taccessions\tresult']
with ThreadPoolExecutor(3) as ex:
    for n, (rid, e, msg) in enumerate(ex.map(one, targets.items())):
        entries[rid] = e; rep.append('\t'.join([rid, targets[rid].get('label', ''), ','.join(targets[rid].get('accs', [])), msg]))
        if (n + 1) % 25 == 0: print('done', n + 1, '/', len(targets), flush=True)
if MISSING:                                                      # merge into the existing file; entries that stay without a map keep the reason found here
    base_ = json.load(open(OUT, encoding='utf-8')) if os.path.exists(OUT) else {'entries': {}}
    base_['entries'].update({k: v for k, v in entries.items() if v.get('loci') or not (base_['entries'].get(k) or {}).get('loci')})
    base_['source'] = (base_.get('source') or '') + ' + NCBI windows for %d entries without a curation-tool item (build_gc_ncbi.py --missing, %s)' % (sum(1 for v in entries.values() if v.get('loci')), time.strftime('%Y-%m-%d'))
    meta = base_
else:
    meta = {'built': time.strftime('%Y-%m-%d'), 'source': 'NCBI (efetch protein ipg + nuccore gb windows)', 'pack': PACK, 'pad_bp': PAD, 'entries': entries}
open(OUT, 'w', encoding='utf-8').write(json.dumps(meta, separators=(',', ':')))
open(os.path.join(REP, 'gc_ncbi_report_%s%s.tsv' % (PACK, '_missing' if MISSING else '')), 'w', encoding='utf-8').write('\n'.join(rep) + '\n')
have = sum(1 for e in entries.values() if e.get('loci'))
print('entries %d, with a genome map %d, without %d; file %d bytes' % (len(entries), have, len(entries) - have, os.path.getsize(OUT)))
for rid, e in entries.items():
    if not e.get('loci'): print('  none:', rid, targets[rid].get('label'), '-', e.get('why', '')[:100])
