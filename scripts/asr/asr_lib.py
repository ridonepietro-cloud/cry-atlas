"""Ancestral sequence reconstruction (ASR) on the trees of the Atlas with PAML (codeml) and GRASP.   Library part: trees, alignment, tool wrappers, result parsing.

Software (installed under atlas/tools/asr, nothing system-wide):
  PAML 4.10.10   https://github.com/abacus-gene/paml      official Windows release (codeml.exe, dat/lg.dat ...)
  GRASP 2.0.0    https://github.com/bodenlab/bnkit        official bnkit.jar (class asr.GRASP; Foley et al. 2022, PLoS Comput Biol), needs Java >= 21: Eclipse Temurin 21 unpacked in tools/asr/jdk
Input = what the page shows: site/data/<pack>/trees.json (IQ-TREE ML tree of a frame x segment set, e.g. struct_A) and the frame alignment (aln_structural.json / aln_sequence.json) cut to the segment and to the columns the tree was made from.
Rooting: the Atlas stores unrooted trees; the reconstruction uses the midpoint root (bifurcating, as GRASP needs); the ancestor of a clade is the node that joins the clade to the rest.
"""
import os, re, sys, json, subprocess, collections, math, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); ATLAS = os.path.dirname(os.path.dirname(HERE))      # project root: packs under site/data
TOOLS = os.environ.get('ASR_TOOLS') or os.path.join(ATLAS, 'tools', 'asr')      # PAML 4.10.10 (Windows build), GRASP bnkit.jar and Temurin JDK 21 unpacked here (set ASR_TOOLS to move them)
CODEML = os.path.join(TOOLS, 'paml', 'paml-4.10.10-win-x86_64', 'bin', 'codeml.exe'); PAML_DAT = os.path.join(TOOLS, 'paml', 'paml-4.10.10-win-x86_64', 'dat')
JAVA = next((os.path.join(TOOLS, 'jdk', d, 'bin', 'java.exe') for d in sorted(os.listdir(os.path.join(TOOLS, 'jdk'))) if d.startswith('jdk')), 'java') if os.path.isdir(os.path.join(TOOLS, 'jdk')) else 'java'
GRASP_JAR = os.path.join(TOOLS, 'bnkit.jar')
AA = 'ACDEFGHIKLMNPQRSTVWY'

# ------------------------------------------------------------------ trees
class Node:
    __slots__ = ('name', 'len', 'kids', 'parent', 'label', 'id')
    def __init__(self, name='', length=0.0, label=''):
        self.name = name; self.len = length; self.kids = []; self.parent = None; self.label = label; self.id = None
def parse_newick(s):
    s = s.strip().rstrip(';'); i = 0; stack = []; cur = Node(); root = cur; n = len(s)
    while i < n:
        ch = s[i]
        if ch == '(': nd = Node(); cur.kids.append(nd); nd.parent = cur; stack.append(cur); cur = nd; i += 1
        elif ch == ',': p = stack[-1]; nd = Node(); p.kids.append(nd); nd.parent = p; cur = nd; i += 1
        elif ch == ')': cur = stack.pop(); i += 1
        elif ch == ':':
            j = i + 1
            while j < n and s[j] not in ',();': j += 1
            cur.len = float(s[i + 1:j]); i = j
        else:
            j = i
            while j < n and s[j] not in ',():;': j += 1
            tok = s[i:j]
            if cur.kids: cur.label = tok
            else: cur.name = tok
            i = j
    return root
def leaves(t):
    out = []; st = [t]
    while st:
        x = st.pop()
        if x.kids: st.extend(reversed(x.kids))
        else: out.append(x)
    return out
def postorder(t):
    out = []; st = [(t, False)]
    while st:
        x, d = st.pop()
        if d: out.append(x)
        else:
            st.append((x, True))
            for c in reversed(x.kids): st.append((c, False))
    return out
def to_newick(t, labels=False, lengths=True):
    def f(x):
        if not x.kids: return x.name + (':%.8f' % x.len if lengths else '')
        return '(' + ','.join(f(c) for c in x.kids) + ')' + (x.label if labels else '') + (':%.8f' % x.len if lengths and x.parent is not None else '')
    return f(t) + ';'
def midpoint_root(t):
    """re-root at the midpoint of the longest leaf-to-leaf path; returns a new bifurcating root (the tree is modified in place)"""
    # distances from every leaf by two passes over the unrooted tree: make adjacency
    nodes = postorder(t); adj = collections.defaultdict(list)
    for x in nodes:
        for c in x.kids: adj[id(x)].append((c, c.len)); adj[id(c)].append((x, c.len))
    byid = {id(x): x for x in nodes}
    def far(src):
        dist = {id(src): 0.0}; prev = {id(src): None}; st = [src]
        while st:
            x = st.pop()
            for y, w in adj[id(x)]:
                if id(y) not in dist: dist[id(y)] = dist[id(x)] + w; prev[id(y)] = x; st.append(y)
        lf = [x for x in nodes if not x.kids]; b = max(lf, key=lambda x: dist[id(x)]); return b, dist, prev
    a = leaves(t)[0]; b, _, _ = far(a); c, dist, prev = far(b); half = dist[id(c)] / 2.0
    # walk from c towards b until the cumulative distance passes half
    x = c
    while prev[id(x)] is not None and dist[id(prev[id(x)])] > half: x = prev[id(x)]
    y = prev[id(x)]; dx = dist[id(x)] - half                                      # the root lies on the edge x - y, dx from x
    # re-root: make a new root on that edge
    ex = x; ey = y; wxy = dist[id(x)] - dist[id(y)]
    # rebuild as a rooted tree by DFS from the new root
    def build(node, came_from):
        n = Node(node.name, 0.0, node.label)
        for nb, w in adj[id(node)]:
            if nb is came_from: continue
            ch = build(nb, node); ch.len = w; n.kids.append(ch); ch.parent = n
        return n
    root = Node(); root.kids = []
    lx = build(ex, ey); lx.len = max(dx, 0.0); ly = build(ey, ex); ly.len = max(wxy - dx, 0.0)
    # the support label of the original edge belongs to the node that was the child side
    lab = ex.label if ex.parent is ey else ey.label
    lx.label = lab if ex.parent is ey else ex.label; ly.label = ey.label if ex.parent is ey else lab
    root.kids = [lx, ly]; lx.parent = root; ly.parent = root
    return root
def leafset(x): return frozenset(l.name for l in leaves(x))
def mrca_node(root, names):
    names = frozenset(names); best = None
    for x in postorder(root):
        ls = leafset(x)
        if names <= ls and (best is None or len(ls) < len(best[1])): best = (x, ls)
    return best[0], best[1]

# ------------------------------------------------------------------ alignment of a tree of the Atlas
def atlas_inputs(pack, key):
    """tree (root node as stored), alignment {leaf: string} cut to the columns of the tree, kept column indices (0-based in the frame), rows"""
    D = os.path.join(ATLAS, 'site', 'data', pack)
    trees = json.load(open(os.path.join(D, 'trees.json'), encoding='utf-8')); meta = json.load(open(os.path.join(D, 'meta.json'), encoding='utf-8'))
    rows = {r['id']: r for r in json.load(open(os.path.join(D, 'rows.json'), encoding='utf-8'))}
    frame, sub = key.split('_'); fr = next(f for f in meta['frames'] if f['id'] == frame)
    a = json.load(open(os.path.join(D, fr['file']), encoding='utf-8')); seqs = a['seqs']; sp = a.get('split', fr['split']); en = fr.get('split2') or fr['width']
    cut = (lambda s: s[:sp]) if sub == 'A' else (lambda s: s[sp:en]) if sub == 'B' else (lambda s: s[:en])
    t = parse_newick(trees['trees'][key]['newick']); ids = [l.name for l in leaves(t)]
    modcols = None
    if sub not in ('A', 'B', 'AB'):                                              # Cry Atlas: a domain tree (DI, DII, DIII) = the columns of that module (the module code of every residue, residues.json) of the frame
        import base64
        mname = meta['module_ids']; res = json.load(open(os.path.join(D, 'residues.json'), encoding='utf-8')); cs = set()
        for rid, g in seqs.items():
            mb = base64.b64decode(res[rid]['A']['mod']); k = 0
            for c, ch in enumerate(g):
                if ch != '-':
                    if mname.get(str(mb[k])) == sub: cs.add(c)
                    k += 1
        modcols = sorted(cs); part = {i: ''.join(seqs[i][c] for c in modcols) for i in ids}
    else:
        part = {i: cut(seqs[i]) for i in ids}
    W = len(next(iter(part.values()))); need = max(1, float(os.environ.get('ASR_MIN_OCC', '0')) * len(ids))        # ASR_MIN_OCC (FliC pack, 2026-11-04: 0.2) = the occupancy of the columns of its trees; default 0 = every column with a residue
    keep = [c for c in range(W) if sum(1 for i in ids if part[i][c] != '-') >= need]
    off = 0 if sub in ('A', 'AB') else sp
    ok = set(AA + '-')                                                           # residues outside the 20 amino acids (X, B, Z, U ...) count as missing, as codeml does; GRASP refuses them
    cols_out = [modcols[c] for c in keep] if modcols is not None else [c + off for c in keep]
    return t, {i: ''.join(part[i][c] if part[i][c] in ok else '-' for c in keep) for i in ids}, cols_out, rows, meta, trees['trees'][key]

# ------------------------------------------------------------------ PAML codeml
PAML_DAT_NAME = {'LG': 'lg', 'WAG': 'wag', 'JTT': 'jones', 'DAYHOFF': 'dayhoff', 'JONES': 'jones'}
def run_paml(root, aln, outdir, model='LG', gamma=True, empirical_freq=True, ncat=4, timeout=int(os.environ.get('ASR_PAML_TIMEOUT', 7200))):      # IRM trees (2,073 leaves): the ancestor step needs more than 2 h (2026-11-02: four runs killed at the limit) -> set ASR_PAML_TIMEOUT
    model = PAML_DAT_NAME.get(model.upper(), model.lower())
    os.makedirs(outdir, exist_ok=True); ids = list(aln); L = len(next(iter(aln.values())))
    with open(os.path.join(outdir, 'aln.phy'), 'w') as f:
        f.write(' %d %d\n' % (len(ids), L))
        for i in ids: f.write('%-12s %s\n' % (i, aln[i]))
    # PAML wants the unrooted tree as a trifurcation or a rooted bifurcating one; the rooted one is given with its branch lengths
    nid = {}; open(os.path.join(outdir, 'tree.nwk'), 'w').write(' %d 1\n%s\n' % (len(ids), to_newick(root, labels=False)))
    ctl = """seqfile = aln.phy
treefile = tree.nwk
outfile = out.txt
noisy = 0
verbose = 0
runmode = 0
seqtype = 2
aaRatefile = %s
model = %d
fix_alpha = 0
alpha = 0.5
Malpha = 0
ncatG = %d
fix_blength = 2
RateAncestor = 1
cleandata = 0
getSE = 0
Small_Diff = 1e-6
method = 0
""" % (model + '.dat', 3 if empirical_freq else 2, ncat if gamma else 1)
    shutil.copyfile(os.path.join(PAML_DAT, model + '.dat'), os.path.join(outdir, model + '.dat'))     # codeml cannot read a path with spaces
    if not gamma: ctl = ctl.replace('fix_alpha = 0', 'fix_alpha = 1').replace('alpha = 0.5', 'alpha = 0')
    open(os.path.join(outdir, 'codeml.ctl'), 'w').write(ctl)
    r = subprocess.run([CODEML, 'codeml.ctl'], cwd=outdir, capture_output=True, text=True, timeout=timeout)
    return r
