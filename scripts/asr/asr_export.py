"""Ancestors of EVERY internal node of the trees of a pack, exported for the Atlas page: site/data/<pack>/asr_<tree key>.json and asr_index.json.

  py asr_export.py <pack> [<tree key> ...]       (all keys of the pack's trees.json when none is given; reconstructs what is not cached: asr_atlas.reconstruct)

File asr_<key>.json (read by the page, lazily):
  { key, pack, tool, model, tree_model, frame, cols [frame columns (0-based) the tree and the ancestors cover], ids [entry ids of the leaves],
    nodes: [ { n: GRASP node label, p: parent node index (-1 root), k: PAML node number, s: support label of the tree, c: leaf indices of the clade, g: GRASP sequence (indel-aware, '-' = deleted),
               a: PAML sequence (every column), q: one digit per column = floor(10 x PAML posterior probability of its residue) capped at 9,
               m: mean PP over the columns GRASP keeps, f: share of those with PP >= 0.9, e: share where PAML = GRASP, x: nearest extant leaf index, y: its identity, r: number of residues } ] }
Rooting: midpoint; a clade of the page whose complement is a clade here has its ancestor at the parent of that complement (the page does that)."""
import os, sys, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import asr_lib as L, asr_parse as P, asr_atlas as A

def export_key(pack, key, indel=None):
    R = A.reconstruct(pack, key, indel=indel); root = R['root']; rows = R['rows']; keep = R['keep']; W = len(keep)
    ids = sorted(R['aln']); idx = {i: k for k, i in enumerate(ids)}; Pm, G = R['paml'], R['grasp']
    pbs = {s: n for n, s in Pm['nodes'].items()}; gbs = {s: n for n, s in G['nodes'].items()}
    post = [x for x in L.postorder(root) if x.kids]; nid = {id(x): k for k, x in enumerate(post)}; out = []
    for x in post:
        ls = L.leafset(x); pn = pbs[ls]; gl = gbs[ls]; pseq = Pm['seq'][pn]; gseq = G['seq'][gl]; pp = Pm['pp'][pn]
        ng = [i for i in range(W) if gseq[i] != '-']; ppv = [pp[i][1] for i in ng]; agree = sum(1 for i in ng if gseq[i] == pseq[i])
        best = (0.0, None)
        for lf in ls:
            a = R['aln'][lf]; m = n = 0
            for i in ng:
                if a[i] != '-': n += 1; m += a[i] == gseq[i]
            if n and m / n > best[0]: best = (m / n, lf)
        out.append({'n': gl, 'p': nid.get(id(x.parent), -1) if x.parent is not None else -1, 'k': pn, 's': x.label or '', 'c': sorted(idx[l] for l in ls), 'g': gseq, 'a': pseq,
                    'q': ''.join(str(min(9, int(p[1] * 10))) for p in pp), 'm': round(sum(ppv) / len(ppv), 3) if ppv else 0, 'f': round(sum(1 for v in ppv if v >= 0.9) / len(ppv), 3) if ppv else 0,
                    'e': round(agree / len(ng), 3) if ng else 0, 'x': idx[best[1]] if best[1] else -1, 'y': round(best[0], 3), 'r': len(ng)})
    meta = R['meta']; fr = next(f for f in meta['frames'] if f['id'] == key.split('_')[0])
    d = {'key': key, 'pack': pack, 'tool': 'PAML 4.10.10 codeml (empirical + F, 4 gamma categories, branch lengths of the tree kept) and GRASP 2.0.0 (indel-aware, indel method ' + R['st'].get('indel', 'BEP') + ', site rates of codeml)', 'indel': R['st'].get('indel', 'BEP'), 'model': R['st']['model'],
         'tree_model': R['tinfo']['model'], 'frame': key.split('_')[0], 'frame_width': fr['width'], 'cols': keep, 'ids': ids, 'nodes': out, 'built': '2026-10-21'}
    path = os.path.join(L.ATLAS, 'site', 'data', pack, 'asr_%s.json' % key)
    json.dump(d, open(path, 'w'), separators=(',', ':'))
    return {'key': key, 'file': 'asr_%s.json' % key, 'nodes': len(out), 'leaves': len(ids), 'columns': W, 'bytes': os.path.getsize(path), 'tree_model': R['tinfo']['model'], 'indel': R['st'].get('indel', 'BEP')}

def main():
    args = sys.argv[1:]; indel = None
    if '--indel' in args: k0 = args.index('--indel'); indel = args[k0 + 1]; del args[k0:k0 + 2]
    sys.argv = [sys.argv[0]] + args
    pack = sys.argv[1]; keys = sys.argv[2:] or list(json.load(open(os.path.join(L.ATLAS, 'site', 'data', pack, 'trees.json'), encoding='utf-8'))['trees'])
    ip = os.path.join(L.ATLAS, 'site', 'data', pack, 'asr_index.json'); idx = json.load(open(ip)) if os.path.exists(ip) else {'pack': pack, 'sets': {}}
    for k in keys:
        try: idx['sets'][k] = export_key(pack, k, indel); print(pack, k, idx['sets'][k], flush=True)
        except Exception as e: print(pack, k, 'FAILED', str(e)[:300], flush=True)
        json.dump(idx, open(ip, 'w'), indent=0)
if __name__ == '__main__':
    main()
