"""Page patch of the Cry Atlas: second choice in the tree window, the NAME-ORDER tree (trees.json ["name_order"], built by 28_name_order_trees.py: a distance tree calculated from the entries in Name order,
order-constrained UPGMA), next to the ML tree.  Needs patch_cry_linked_tree.py applied first (the ML tree with lines to rows in another order); applied IN PLACE to site/index.html, idempotent.
What it changes:
  * Tree window selector: 'ML tree' (default) and 'Name-order tree'; the frame (structure / sequence) picks the matching tree.
  * Choosing the name-order tree puts the rows in Name order, where tips and rows are level (no lines); with another order its tips are joined to the rows by lines like the ML tree.
    Choosing the ML tree while the rows are listed by Name keeps the Name order (tips joined by lines), as before; in other orders it puts the rows in tree order. Order = ML tree switches the window back to the ML tree.
  * Order = Name uses the order stored with the trees (the page's own case-insensitive sort) so that tips and rows cannot disagree between browser locales (old sort used when it is not stored).
  * The tree is drawn exactly as stored (no re-rooting, no ladderizing); there are no ancestors for it (clicking a clade picks the entries only); Newick download is named *_name_order.
Run: py -I patch_cry_name_tree.py [path of index.html]"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); CRY = os.path.dirname(HERE)
P = sys.argv[1] if len(sys.argv) > 1 else os.path.join(CRY, "site", "index.html")
t = open(P, encoding="utf-8", newline="").read()
MARK = "/* name-order tree (patch_cry_name_tree.py) */"
if MARK in t:
    print("already patched:", P); sys.exit(0)
assert "/* linked tree window (patch_cry_linked_tree.py) */" in t, "apply patch_cry_linked_tree.py first"


def rep(old, new):
    global t
    assert t.count(old) == 1, "expected exactly one match (%d) for: %s" % (t.count(old), old[:100])
    t = t.replace(old, new)


# 1. Order = Name: the order stored with the trees (same as the page's own sort), then the label as tie-break
rep('if (S.order === "name") rows.sort(function (a, b) { return a.label.localeCompare(b.label); });',
    'if (S.order === "name") { var nrk = nameRank(); rows.sort(function (a, b) { return nrk ? ((nrk[a.id] == null ? 1e9 : nrk[a.id]) - (nrk[b.id] == null ? 1e9 : nrk[b.id])) || a.label.localeCompare(b.label, "en") : a.label.localeCompare(b.label); }); }')

# 2. helpers (before tdKey); the order in which the tree on show stands level with its rows
helpers = MARK + r'''
function nameRank() {                                          // rank of every entry in the Name order stored with the trees (null: not stored)
  var no = S.trees && S.trees.name_order && S.trees.name_order.order; if (!no) return null;
  if (S._nameRankSrc !== no) { var m = {}; no.forEach(function (id, i) { m[id] = i; }); S._nameRank = m; S._nameRankSrc = no; }
  return S._nameRank;
}
function nameParse(str) {                                      // rooted Newick -> {name, len, sup, kids}, child order kept (the leaf order IS the Name order)
  var i = 0;
  function token() { var s = i; while (i < str.length && ",():;".indexOf(str[i]) < 0) i++; return str.slice(s, i); }
  function sub() {
    var o = { name: "", len: 0, sup: null, kids: [] };
    if (str[i] === "(") { i++; for (;;) { o.kids.push(sub()); if (str[i] === ",") { i++; continue; } if (str[i] === ")") { i++; break; } throw new Error("Newick: unexpected character at " + i); } token(); }
    else o.name = token();
    if (str[i] === ":") { i++; o.len = parseFloat(token()); }
    return o;
  }
  var r = sub(); r.len = 0; return r;
}
function tdWantOrder() { return TD.set === "N" ? "name" : "ml:" + TD.set; }     // row order in which the tree on show stands level with its rows
function tdFollowOrder() {                                     // the Order selector changed: tree order -> the ML tree again
  if (!TD.on || TD.set !== "N" || S.order.indexOf("ml:") !== 0) return;
  TD.set = "A"; TD.built = null; tdSave(); tdSyncUi();
}
function tdKey() {'''
rep("function tdKey() {", helpers)
rep('function tdLinked() { return S.order !== "ml:" + TD.set; }', "function tdLinked() { return S.order !== tdWantOrder(); }")

# 3. the tree as stored
rep("  var key = tdKey(); if (!S.trees || !S.trees.trees[key]) return null;\n  if (!TD.root[key])",
    '  if (TD.set === "N") { var nt = S.trees && S.trees.name_order && S.trees.name_order[S.frame]; if (!nt) return null; TD.nroot = TD.nroot || {};\n'
    '    if (!TD.nroot[S.frame] || TD.nroot[S.frame].nw !== nt.newick) TD.nroot[S.frame] = { nw: nt.newick, root: nameParse(nt.newick) };\n    return TD.nroot[S.frame].root; }\n'
    "  var key = tdKey(); if (!S.trees || !S.trees.trees[key]) return null;\n  if (!TD.root[key])")

# 4. footer, download, ancestors
rep("var ti = S.trees.trees[tdKey()];", 'var ti = TD.set === "N" ? null : S.trees.trees[tdKey()];')
rep('"midpoint root", 6, h - BOT + 17);', '(TD.set === "N" ? "name-order tree (UPGMA, not ML)" : "midpoint root"), 6, h - BOT + 17);')
rep('dlSaveTree(root, tdKey() + "_midpoint", withNames);', 'dlSaveTree(root, tdKey() + (TD.set === "N" ? "_name_order" : "_midpoint"), withNames);')
rep('"Newick of the tree on show: the entries on show, midpoint-rooted; branch lengths in substitutions per site"',
    '"Newick of the tree on show: the entries on show, " + (TD.set === "N" ? "rooted as drawn (leaves in name order)" : "midpoint-rooted") + "; branch lengths in substitutions per site"')
rep("if (!e.altKey) anOnPick(ids, mlOrderKey(TD.set));", 'if (!e.altKey && TD.set !== "N") anOnPick(ids, mlOrderKey(TD.set));')

# 5. ordering and switching
rep('function tdSyncOrder() { S.order = "ml:" + TD.set;', "function tdSyncOrder() { S.order = tdWantOrder();")
rep('if (TD.on && S.order !== "ml:" + TD.set && !tdKeepOrder()) tdSyncOrder(); else { sizeScroller(); schedule(); }', 'if (TD.on && S.order !== tdWantOrder() && !tdKeepOrder()) tdSyncOrder(); else { sizeScroller(); schedule(); }')
rep('if (TD.on && S.order !== "ml:" + TD.set && !tdKeepOrder()) { S.order = "ml:" + TD.set; $("orderSel").value = S.order; }', 'if (TD.on && S.order !== tdWantOrder() && !tdKeepOrder()) { S.order = tdWantOrder(); $("orderSel").value = S.order; }')
rep('if (["A", "B", "AB"].indexOf(TD.set) < 0 || S.meta.homo) TD.set = S.meta.homo ? "A" : "AB";',
    'if (!(S.meta.homo && TD.set === "N" && S.trees && S.trees.name_order) && (["A", "B", "AB"].indexOf(TD.set) < 0 || S.meta.homo)) TD.set = S.meta.homo ? "A" : "AB";\n'
    '  if (S.trees && S.trees.name_order && S.meta.homo && $("tdSet") && !$("tdSet").querySelector(\'option[value="N"]\')) {\n'
    '    var tsel = $("tdSet"); tsel.options[0].textContent = "ML tree"; tsel.style.maxWidth = "130px"; tsel.title = "ML tree: maximum-likelihood tree of the frame on screen. Name-order tree: distance tree (UPGMA) calculated from the entries in name order, only name neighbours can be joined; its tips are level with the rows when Order = Name";\n'
    '    var nop = document.createElement("option"); nop.value = "N"; nop.textContent = "Name-order tree"; tsel.appendChild(nop);\n  }')
rep('$("tdSet").addEventListener("change", function (e) { TD.set = e.target.value; tdSave(); TD.built = null; tdSyncOrder(); });',
    '$("tdSet").addEventListener("change", function (e) { TD.set = e.target.value; tdSave(); TD.built = null; if (TD.set !== "N" && tdKeepOrder()) { tdSyncUi(); sizeScroller(); schedule(); } else tdSyncOrder(); });')
rep('$("orderSel").addEventListener("change", function (e) { S.order = e.target.value; store("order", S.order); rebuild(); });',
    '$("orderSel").addEventListener("change", function (e) { S.order = e.target.value; store("order", S.order); tdFollowOrder(); rebuild(); });')

open(P, "w", encoding="utf-8", newline="").write(t)
print("patched", P, len(t), "bytes")
