"""Ancestral sequence reconstruction of nodes of the Atlas trees with PAML (codeml) and GRASP.

  py asr_atlas.py <pack> <tree key> [--default] [--class NAME] [--clade id1,id2,...] [--name NAME] [--model LG] [--force]

  pack      fit | git262 | gldlm | gitx | mscl | mscs | piezo          tree key   struct_A | struct_B | struct_AB | seq_A | seq_B | seq_AB   (the sets of site/data/<pack>/trees.json)
  --default the nodes of interest of the pack: for every class the largest clade with >= 90 % of that class (>= 5 entries), the tree root of the midpoint rooting, and the two clades below its larger child
  --class   the smallest clade that contains every entry of that class        --clade  the smallest clade that contains these entries (row ids)

The whole tree is reconstructed once per pack / tree key / model (codeml: empirical matrix + F, 4 gamma categories, branch lengths of the IQ-TREE tree kept; GRASP: same matrix, site rates from codeml, indels by bi-directional edge parsimony)
and cached in out/<pack>_<key>/ ; every node asked for is then written to out/<pack>_<key>/nodes/ :
   <name>.fa           ancestors in the columns of the Atlas frame (gaps where the tree has no column):  paml (every column), grasp (indel-aware), and the extant entry most similar to the grasp ancestor
   <name>_columns.tsv  per frame column: GRASP residue, PAML residue and posterior probability, whether they agree
   nodes.tsv / ancestors.json  one line / record per node with composition, support, mean posterior probability, agreement of the two programs.
Software: asr_lib.py (PAML 4.10.10, GRASP 2.0.0 in atlas/tools/asr)."""
import os, sys, json, time, collections, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import asr_lib as L, asr_parse as P

OUT = os.path.join(L.HERE, 'out')
def reconstruct(pack, key, model='LG', force=False, threads=4, indel=None):
    d = os.path.join(OUT, '%s_%s' % (pack, key)); os.makedirs(d, exist_ok=True)
    t, aln, keep, rows, meta, tinfo = L.atlas_inputs(pack, key); root = L.midpoint_root(t)
    st = {'pack': pack, 'key': key, 'model': model, 'leaves': len(aln), 'columns': len(keep), 'tree_model': tinfo['model']}
    pd, gd = os.path.join(d, 'paml'), os.path.join(d, 'grasp')
    if force or not os.path.exists(os.path.join(pd, 'rst')):
        t0 = time.time(); r = L.run_paml(root, aln, pd, model=model, gamma=True, empirical_freq=True)
        if not os.path.exists(os.path.join(pd, 'rst')): raise RuntimeError('codeml failed: ' + r.stdout[-400:])
        st['paml_s'] = round(time.time() - t0)
    rates = P.read_paml_rates(os.path.join(pd, 'rates'))
    if force or not os.path.exists(os.path.join(gd, 'out', 'aln_ancestors.fa')):
        t0 = time.time(); r = P.run_grasp(root, aln, gd, model=model, threads=int(os.environ.get('ASR_GRASP_THREADS', threads)), mem=os.environ.get('ASR_GRASP_MEM', '8g'), rates=rates if len(rates) == len(keep) else None, indel=indel, chunks=int(os.environ.get('ASR_GRASP_CHUNKS', 1)))
        if not os.path.exists(os.path.join(gd, 'out', 'aln_ancestors.fa')): raise RuntimeError('GRASP failed: ' + (r.stdout + r.stderr)[-600:])
        st['grasp_s'] = round(time.time() - t0)
        open(os.path.join(gd, 'indel.txt'), 'w').write((indel or 'BEP') + (', columns in %s blocks' % os.environ['ASR_GRASP_CHUNKS'] if int(os.environ.get('ASR_GRASP_CHUNKS', 1)) > 1 else ''))
    ip = os.path.join(gd, 'indel.txt'); st['indel'] = open(ip).read().strip() if os.path.exists(ip) else 'BEP'
    paml = P.parse_rst(os.path.join(pd, 'rst')); grasp = P.parse_grasp(gd)
    return dict(root=root, aln=aln, keep=keep, rows=rows, meta=meta, tinfo=tinfo, paml=paml, grasp=grasp, st=st, dir=d, key=key, pack=pack)

def pure_clades(root, rows, cls, min_purity=0.9, min_size=5):
    """maximal clades in which >= min_purity of the entries are of class cls (and that contain >= min_size of them)"""
    best = []
    def walk(x):
        ls = leafset(x); n = sum(1 for i in ls if rows[i]['class'] == cls)
        if n >= min_size and n / len(ls) >= min_purity: best.append((n, x)); return
        for c in x.kids: walk(c)
    walk(root); best.sort(key=lambda a: -a[0]); return best
leafset = L.leafset

def resolve(R, spec):
    """node spec -> (node, name)"""
    root, rows = R['root'], R['rows']
    if spec.get('clade'): x, ls = L.mrca_node(root, spec['clade']); return x, spec.get('name') or 'clade_%d' % len(ls)
    if spec.get('class'):
        ids = [l.name for l in L.leaves(root) if rows[l.name]['class'] == spec['class']]; x, ls = L.mrca_node(root, ids); return x, spec.get('name') or 'MRCA_%s' % spec['class']
    raise ValueError(spec)

def default_nodes(R):
    root, rows = R['root'], R['rows']; out = []; seen = set()
    def add(x, name):
        ls = leafset(x)
        if ls in seen or len(ls) < 5: return
        seen.add(ls); out.append((x, name))
    add(root, 'root')
    classes = collections.Counter(rows[l.name]['class'] for l in L.leaves(root))
    for c, n in classes.most_common():
        if n < 8: continue
        pc = pure_clades(root, rows, c)
        for k, (m, x) in enumerate(pc[:2]): add(x, 'clade_%s%s' % (c.replace(' ', '_').replace('(', '').replace(')', ''), '' if k == 0 else '_%d' % (k + 1)))
    big = max(root.kids, key=lambda k: len(leafset(k)))
    add(big, 'main_clade')
    for j, k in enumerate(sorted(big.kids, key=lambda k: -len(leafset(k)))): add(k, 'main_%d' % (j + 1))
    return out

def extract(R, x, name):
    """sequences and statistics of one node"""
    ls = leafset(x); P_, G = R['paml'], R['grasp']; rows = R['rows']
    pn = next((n for n, s in P_['nodes'].items() if s == ls), None); gl = next((n for n, s in G['nodes'].items() if s == ls), None)
    if pn is None or gl is None: raise RuntimeError('node %s not found in the output of %s' % (name, 'PAML' if pn is None else 'GRASP'))
    pseq = P_['seq'][pn]; gseq = G['seq'][gl]; pp = P_['pp'][pn]; keep = R['keep']; W = len(keep)
    assert len(pseq) == W and len(gseq) == W and len(pp) == W, (len(pseq), len(gseq), len(pp), W)
    nongap = [i for i in range(W) if gseq[i] != '-']
    agree = [i for i in nongap if gseq[i] == pseq[i]]
    ppv = [pp[i][1] for i in nongap]
    # nearest extant by identity over the columns the GRASP ancestor occupies
    best = (0, None)
    for lf in ls:
        a = R['aln'][lf]; m = n = 0
        for i in nongap:
            if a[i] != '-': n += 1; m += a[i] == gseq[i]
        idn = m / n if n else 0
        if idn > best[0]: best = (idn, lf)
    comp = collections.Counter(rows[i]['class'] for i in ls)
    rec = {'name': name, 'pack': R['pack'], 'key': R['key'], 'n_leaves': len(ls), 'classes': dict(comp), 'support': x.label, 'paml_node': pn, 'grasp_node': gl,
           'residues': len(nongap), 'columns': W, 'mean_pp': round(sum(ppv) / len(ppv), 3) if ppv else None, 'frac_pp90': round(sum(1 for v in ppv if v >= 0.9) / len(ppv), 3) if ppv else None,
           'frac_pp50': round(sum(1 for v in ppv if v >= 0.5) / len(ppv), 3) if ppv else None,
           'agree_paml_grasp': round(len(agree) / len(nongap), 3) if nongap else None, 'nearest_extant': best[1], 'nearest_identity': round(best[0], 3)}
    def full(s):                                                                       # onto the columns of the Atlas frame
        f = ['-'] * R['frame_width']
        for k, c in enumerate(keep): f[c] = s[k]
        return ''.join(f)
    rec['paml'] = pseq; rec['grasp'] = gseq; rec['paml_pp'] = [round(v[1], 3) for v in pp]; rec['kept_columns'] = keep
    return rec, pseq, gseq, pp, best[1]

def write_node(R, rec, pseq, gseq, pp, nearest):
    nd = os.path.join(R['dir'], 'nodes'); os.makedirs(nd, exist_ok=True); name = rec['name']; keep = R['keep']; aln = R['aln']
    def full(s):
        f = ['-'] * R['frame_width']
        for k, c in enumerate(keep): f[c] = s[k]
        return ''.join(f)
    with open(os.path.join(nd, name + '.fa'), 'w') as f:
        f.write('>%s_grasp node %s (%d entries)\n%s\n' % (name, rec['grasp_node'], rec['n_leaves'], full(gseq)))
        f.write('>%s_paml node %d mean PP %s\n%s\n' % (name, rec['paml_node'], rec['mean_pp'], full(pseq)))
        if nearest: f.write('>%s_nearest_extant %s identity %.3f\n%s\n' % (name, nearest, rec['nearest_identity'], full(aln[nearest])))
    with open(os.path.join(nd, name + '_columns.tsv'), 'w') as f:
        f.write('frame_column\tgrasp\tpaml\tpaml_pp\tagree\n')
        for k, c in enumerate(keep): f.write('%d\t%s\t%s\t%.3f\t%s\n' % (c + 1, gseq[k], pseq[k], pp[k][1], '' if gseq[k] == '-' else int(gseq[k] == pseq[k])))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('pack'); ap.add_argument('key'); ap.add_argument('--default', action='store_true'); ap.add_argument('--class', dest='cls'); ap.add_argument('--clade'); ap.add_argument('--name')
    ap.add_argument('--model', default='LG'); ap.add_argument('--force', action='store_true'); ap.add_argument('--threads', type=int, default=4); a = ap.parse_args()
    R = reconstruct(a.pack, a.key, a.model, a.force, a.threads)
    meta = R['meta']; fr = next(f for f in meta['frames'] if f['id'] == a.key.split('_')[0]); R['frame_width'] = fr['width']
    print('reconstructed', R['st'])
    nodes = []
    if a.default: nodes += default_nodes(R)
    if a.cls: nodes.append(resolve(R, {'class': a.cls, 'name': a.name}))
    if a.clade: nodes.append(resolve(R, {'clade': a.clade.split(','), 'name': a.name}))
    recs = []
    for x, name in nodes:
        rec, ps, gs, pp, near = extract(R, x, name); write_node(R, rec, ps, gs, pp, near); recs.append(rec)
        print('%-26s %4d entries  support %-5s residues %4d  mean PP %.3f  PP>=0.9 %.0f %%  PAML = GRASP %.0f %%  nearest %s %.0f %%  %s' % (name, rec['n_leaves'], rec['support'], rec['residues'], rec['mean_pp'], 100 * rec['frac_pp90'], 100 * rec['agree_paml_grasp'], rec['nearest_extant'], 100 * rec['nearest_identity'], rec['classes']))
    p = os.path.join(R['dir'], 'ancestors.json'); old = json.load(open(p)) if os.path.exists(p) else {}
    for r in recs: old[r['name']] = r
    json.dump({'info': R['st'], 'nodes': old}, open(p, 'w'))
    with open(os.path.join(R['dir'], 'nodes.tsv'), 'w') as f:
        f.write('name\tentries\tsupport\tresidues\tmean_pp\tfrac_pp90\tpaml_eq_grasp\tnearest_extant\tnearest_identity\tclasses\n')
        for r in old.values(): f.write('%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n' % (r['name'], r['n_leaves'], r['support'], r['residues'], r['mean_pp'], r['frac_pp90'], r['agree_paml_grasp'], r['nearest_extant'], r['nearest_identity'], json.dumps(r['classes'])))
if __name__ == '__main__':
    main()
