"""Readers of the outputs of codeml (rst, rates) and GRASP (aln_ancestors.fa / .nwk), and the GRASP runner.  Part of the ASR pipeline (asr_lib.py, asr_atlas.py)."""
import os, re, subprocess
from asr_lib import JAVA, GRASP_JAR, to_newick, parse_newick, postorder, leafset

def parse_rst(path):
    """codeml rst -> {'ntips', 'nodes': {node number: frozenset(leaf names)}, 'seq': {node number: sequence}, 'pp': {node number: [(aa, prob) per site]}, 'first', 'last'}"""
    txt = open(path, errors='replace').read()
    m = re.search(r"tree with node labels for Rod Page's TreeView\s*\n(.*?;)", txt, re.S)
    toks = re.findall(r"\(|\)|,|;|\d+_[^,():;\s]+|\d+", m.group(1))
    stack = []; cur = []
    for t in toks:
        if t == '(': stack.append(cur); cur = []
        elif t == ')': kids = cur; cur = stack.pop(); cur.append(('int', kids))
        elif t in ',;': pass
        elif '_' in t and t.split('_')[0].isdigit(): cur.append(('tip', t.split('_', 1)[1]))
        else: cur[-1] = ('int', cur[-1][1], int(t))                                # a number after ")" labels the internal node
    root = cur[-1]; sets = {}
    def ls(x):
        if x[0] == 'tip': return frozenset([x[1]])
        r = frozenset().union(*[ls(c) for c in x[1]])
        if len(x) > 2: sets[x[2]] = r
        return r
    allt = ls(root)
    seqs = {}
    sec = txt.split('List of extant and reconstructed sequences', 1)[1]
    for ln in sec.splitlines():
        if ln.startswith('node #'):
            parts = ln.split(); seqs[int(parts[1].lstrip('#'))] = ''.join(parts[2:])
    mm = re.search(r'Nodes (\d+) to (\d+) are ancestral', txt); first, last = int(mm.group(1)), int(mm.group(2))
    body = txt.split('Prob of best state at each node, listed by site', 1)[1].split('Summary of changes along branches', 1)[0]
    pp = {n: [] for n in range(first, last + 1)}
    for ln in body.splitlines():
        m2 = re.match(r'^\s*(\d+)\s+(\d+)\s+(\S+)\s*:\s*(.*)$', ln)
        if not m2: continue
        for k, (aa, pr) in enumerate(re.findall(r'(\S)\(([0-9.]+)\)', m2.group(4))): pp[first + k].append((aa, float(pr)))
    return {'ntips': len(allt), 'nodes': sets, 'seq': seqs, 'pp': pp, 'first': first, 'last': last}

def read_paml_rates(path):
    """mean posterior rate per site from the codeml file 'rates'"""
    out = []
    for ln in open(path):
        m = re.match(r'^\s*\d+\s+\d+\s+\S+\s+([0-9.]+)\s+\d+\s*$', ln)
        if m: out.append(float(m.group(1)))
    return out

def run_grasp(root, aln, outdir, model='LG', threads=4, rates=None, mem='8g', timeout=7200, indel=None, chunks=1):
    """GRASP on the alignment; chunks > 1 splits the columns into blocks that are run one after the other on the same tree and joined by node label (position-specific indel methods such as PSP
    and the column-wise character inference treat the columns independently; needed for trees of several thousand columns, where the character inference of GRASP does not fit in 8 GB)"""
    os.makedirs(outdir, exist_ok=True)
    if chunks > 1:                                                           # Cry Atlas change: INTERLEAVED blocks (columns c, c + chunks, c + 2 chunks ...) instead of contiguous ones: a contiguous block can lie in a region that a whole clade lacks
        W = len(next(iter(aln.values()))); blocks = []; res = None            # (the protoxin half of Cry2 / Cry3 ...), and GRASP then finds no path for the ancestors of that clade; valid for column-independent indel methods (PSP)
        for c in range(chunks):
            sub = {k: v[c::chunks] for k, v in aln.items()}
            if not len(next(iter(sub.values()))): break
            res = run_grasp(root, sub, os.path.join(outdir, 'chunk_%d' % c), model=model, threads=threads, rates=(rates[c::chunks] if rates else None), mem=mem, timeout=timeout, indel=indel, chunks=1)
            fp = os.path.join(outdir, 'chunk_%d' % c, 'out', 'aln_ancestors.fa')
            if not os.path.exists(fp): return res
            seqs = {}; cur = None
            for ln in open(fp):
                ln = ln.rstrip()
                if ln.startswith('>'): cur = ln[1:].strip(); seqs.setdefault(cur, [])
                elif cur is not None: seqs[cur].append(ln)
            blocks.append({k: ''.join(v) for k, v in seqs.items()})
            if c == 0: nwk0 = open(os.path.join(outdir, 'chunk_0', 'out', 'aln_ancestors.nwk')).read()
            else: assert open(os.path.join(outdir, 'chunk_%d' % c, 'out', 'aln_ancestors.nwk')).read() == nwk0, 'node labels differ between blocks'
        os.makedirs(os.path.join(outdir, 'out'), exist_ok=True)
        with open(os.path.join(outdir, 'out', 'aln_ancestors.fa'), 'w') as f:
            for k in blocks[0]:
                parts = [b[k] for b in blocks]
                assert all(len(parts[i]) == len(range(i, W, chunks)) for i in range(len(parts))), 'block length differs for ' + k
                f.write('>' + k + chr(10) + ''.join(parts[j % chunks][j // chunks] for j in range(W)) + chr(10))
        open(os.path.join(outdir, 'out', 'aln_ancestors.nwk'), 'w').write(nwk0)
        return res
    with open(os.path.join(outdir, 'aln.fa'), 'w') as f:
        for i, sq in aln.items(): f.write('>%s\n%s\n' % (i, sq))
    open(os.path.join(outdir, 'tree.nwk'), 'w').write(to_newick(root, labels=False) + '\n')
    cmd = [JAVA, '-Xmx' + mem, '-jar', GRASP_JAR, '--aln', 'aln.fa', '--nwk', 'tree.nwk', '--output-folder', 'out', '--substitution-model', model, '--save-as', 'FASTA', 'TREE', '--threads', str(threads), '--verbose', '--time']
    if indel and indel != 'BEP': cmd += ['--indel-method', indel]                      # BEP (default) did not finish in hours on the A+B trees of the GIT pack (2,556 columns); PSP takes < 1 min
    if rates:
        with open(os.path.join(outdir, 'rates.tsv'), 'w') as f:
            f.write('Site\tRate\n')
            for k, v in enumerate(rates, 1): f.write('%d\t%.5f\n' % (k, v))
        cmd += ['--rates-file', 'rates.tsv']
    return subprocess.run(cmd, cwd=outdir, capture_output=True, text=True, timeout=timeout)

def parse_grasp(outdir):
    """-> {'nodes': {N label: frozenset(leaf names)}, 'seq': {N label: gapped aligned sequence}}"""
    fa = {}; name = None
    for ln in open(os.path.join(outdir, 'out', 'aln_ancestors.fa')):
        ln = ln.rstrip()
        if ln.startswith('>'): name = ln[1:].strip(); fa[name] = []
        elif name: fa[name].append(ln)
    seq = {k: ''.join(v) for k, v in fa.items()}
    t = parse_newick(open(os.path.join(outdir, 'out', 'aln_ancestors.nwk')).read()); sets = {}
    for x in postorder(t):
        if x.kids and x.label: sets[x.label] = leafset(x)
    return {'nodes': sets, 'seq': seq}
