"""Page patch of the Cry Atlas: the tree window beside the alignment also works when the rows are NOT in the order of the tree (e.g. Order = Name).
The ML tree keeps its own (ladderized, midpoint-rooted) leaf order; every tip is joined by a line to its row, so the crossings show how far the order of the rows differs from the tree.
Applied IN PLACE to site/index.html after patch_cry_page.py and the downloads patch (idempotent: a marker comment stops a second run).
What it changes (tree window = 'dock', functions tdBuild / drawTreeDock / tdPointer / tdToggle / initTreeDock):
  * rows in another order than ML tree (Name, tree order of the guide tree, curation orders): instead of the note 'The rows are not in the order of this tree' the tree is drawn at its own height
    (tip k at the k-th row position, same row height and scrolling as the rows) and a line joins each tip to its row, coloured by family; the entries picked / on screen / hovered are drawn bold;
    the button 'Rows by tree' in the control bar of the window puts the rows in tree order (as before). Hovering a clade or a tip highlights the lines and rows of its entries; clicking works as before.
  * opening the window no longer re-sorts rows that are listed by Name (or by a curation order); other orders are still switched to the ML tree order when the window opens.
Run: py -I patch_cry_linked_tree.py [path of index.html]"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); CRY = os.path.dirname(HERE)
P = sys.argv[1] if len(sys.argv) > 1 else os.path.join(CRY, "site", "index.html")
t = open(P, encoding="utf-8", newline="").read()
MARK = "/* linked tree window (patch_cry_linked_tree.py) */"
if MARK in t:
    print("already patched:", P); sys.exit(0)


def rep(old, new):
    global t
    assert t.count(old) == 1, "expected exactly one match (%d) for: %s" % (t.count(old), old[:100])
    t = t.replace(old, new)


helpers = MARK + r'''
function tdLinked() { return S.order !== "ml:" + TD.set; }            // rows are not in the order of the tree: tips are joined to their rows by lines
function tdKeepOrder() { return S.order === "name" || S.order.indexOf("cu:") === 0; }   // orders the tree window must not replace when it opens
function tdLinks(ctx, B, V, x0, x1, rh, st, h, hr, hov) {            // one line per tip: (x0, tip position) -> (x1, row position); bold for picked, current and hovered entries
  var groups = {}, bold = [], vTop = TOP - 2, vBot = h - BOT + 2;
  B.leaves.forEach(function (l) {
    var ty = TOP + (l.y + 0.5) * rh - st, ry = TOP + (l.ri + 0.5) * rh - st;
    if ((ty < vTop && ry < vTop) || (ty > vBot && ry > vBot)) return;
    var r = V.rows[l.ri], cur = S.cur && S.cur.id === r.id, pk = S.picked.has(r.id), hot = hr === l.ri || (hov && l.y >= hov.lo && l.y <= hov.hi);
    if (cur || pk || hot) bold.push({ ty: ty, ry: ry, cur: cur, hot: hot }); else { var c = tdTipColor(r); (groups[c] = groups[c] || []).push(ty, ry); }
  });
  ctx.lineWidth = 1; ctx.globalAlpha = rh >= 6 ? 0.6 : 0.4;
  Object.keys(groups).forEach(function (c) { ctx.strokeStyle = c; ctx.beginPath(); var a = groups[c]; for (var i = 0; i < a.length; i += 2) { ctx.moveTo(x0, a[i]); ctx.lineTo(x1, a[i + 1]); } ctx.stroke(); });
  ctx.globalAlpha = 1;
  bold.forEach(function (e) { ctx.strokeStyle = e.cur || e.hot ? T.accent : T["cv-ink"]; ctx.lineWidth = e.cur ? 2.2 : 1.7; ctx.beginPath(); ctx.moveTo(x0, e.ty); ctx.lineTo(x1, e.ry); ctx.stroke(); });
}
function tdKey() {'''
rep("function tdKey() {", helpers)

# layout: tip position of the tree (linked) or row position (aligned); leaf.ri = row index of the leaf's entry
rep("(function y(o) { if (!o.kids.length) { o.y = V.rowIdx[o.name]; out.leaves.push(o); return; } o.kids.forEach(y); o.y = (o.kids[0].y + o.kids[o.kids.length - 1].y) / 2; })(root);",
    "out.link = tdLinked(); out.byRow = []; (function y(o) { if (!o.kids.length) { o.ri = V.rowIdx[o.name]; o.y = out.link ? out.leaves.length : o.ri; out.byRow[o.ri] = o; out.leaves.push(o); return; } o.kids.forEach(y); o.y = (o.kids[0].y + o.kids[o.kids.length - 1].y) / 2; })(root);")

# drawing
rep('  if (S.order !== "ml:" + TD.set) { tdNote(ctx, w, ["The rows are not in the order", "of this tree."], "Order rows by this tree"); return; }\n',
    '  var rb0 = $("tdRowsBy"); if (rb0) rb0.hidden = !tdLinked();                  // button of the control bar: only while the rows are in another order\n')
rep("TD.built.key !== tdKey() || TD.built.lens !== TD.lens) TD.built = tdBuild();", "TD.built.key !== tdKey() || TD.built.lens !== TD.lens || TD.built.link !== tdLinked()) TD.built = tdBuild();")
rep("var padL = 8, padR = 10, treeW = Math.max(40, w - padL - padR - 8), X = function (o) { return padL + o.d / B.maxD * treeW; }, Y = function (o) { return TOP + (o.y + 0.5) * rh - st; };",
    "var padL = 8, padR = 10, treeW = Math.max(40, (w - padL - padR - 8) * (B.link ? 0.58 : 1)), X = function (o) { return padL + o.d / B.maxD * treeW; }, Y = function (o) { return TOP + (o.y + 0.5) * rh - st; };\n"
    "  var x0 = padL + treeW + 8, bx = B.link ? x0 - 2 : 0; TD.linkX = B.link ? x0 : null;")
rep('if (S.picked.has(rr.id)) { ctx.fillStyle = T["cv-pick"]; ctx.fillRect(0, y0, w, rh); }', 'if (S.picked.has(rr.id)) { ctx.fillStyle = T["cv-pick"]; ctx.fillRect(bx, y0, w - bx, rh); }')
rep('if (S.cur && S.cur.id === rr.id) { ctx.fillStyle = T["cv-hl"]; ctx.fillRect(0, y0, w, rh); ctx.fillStyle = T.accent; ctx.fillRect(w - 3, y0, 3, rh); }',
    'if (S.cur && S.cur.id === rr.id) { ctx.fillStyle = T["cv-hl"]; ctx.fillRect(bx, y0, w - bx, rh); ctx.fillStyle = T.accent; ctx.fillRect(w - 3, y0, 3, rh); }')
rep('if (hov) { ctx.fillStyle = T["cv-hl"]; ctx.fillRect(0, rowY(hov.lo), w, (hov.hi - hov.lo + 1) * rh); }',
    'if (hov) { ctx.fillStyle = T["cv-hl"]; if (B.link) { ctx.fillRect(0, rowY(hov.lo), x0 - 2, (hov.hi - hov.lo + 1) * rh); B.leaves.forEach(function (l) { if (l.y >= hov.lo && l.y <= hov.hi) ctx.fillRect(bx, rowY(l.ri), w - bx, rh); }); } else ctx.fillRect(0, rowY(hov.lo), w, (hov.hi - hov.lo + 1) * rh); }')
rep('if (hr != null) { ctx.fillStyle = T["cv-hl"]; ctx.fillRect(0, rowY(hr), w, rh); }',
    'if (hr != null) { ctx.fillStyle = T["cv-hl"]; if (B.link) { ctx.fillRect(bx, rowY(hr), w - bx, rh); if (B.byRow[hr]) ctx.fillRect(0, rowY(B.byRow[hr].y), x0 - 2, rh); } else ctx.fillRect(0, rowY(hr), w, rh); }')
rep("  B.leaves.forEach(function (l) {\n    if (l.y < r0 - 1 || l.y > r1 + 1) return;\n    var r = V.rows[l.y], x = X(l)",
    "  if (B.link) tdLinks(ctx, B, V, x0, w - 3, rh, st, h, hr, hov);\n  B.leaves.forEach(function (l) {\n    if (l.y < r0 - 1 || l.y > r1 + 1) return;\n    var r = V.rows[l.ri], x = X(l)")
rep("ctx.lineTo(w - padR + 2, y); ctx.stroke(); ctx.setLineDash([]); }", "ctx.lineTo(B.link ? x0 - 2 : w - padR + 2, y); ctx.stroke(); ctx.setLineDash([]); }")
rep("Math.min(pk || S.rowHot === l.y ? 4 : 3.4,", "Math.min(pk || S.rowHot === l.ri ? 4 : 3.4,")

# pointer: left of the lines the pointer is over a tip position, right of them over a row
rep("var ri = Math.floor((y - TOP + st) / S.rh); return ri >= 0 && ri < V.nr ? { ri: ri } : null;",
    "var ri = Math.floor((y - TOP + st) / S.rh); if (TD.linkX != null && TD.built && TD.built.link && x < TD.linkX) { var lf = TD.built.leaves[ri]; return lf ? { ri: lf.ri } : null; } return ri >= 0 && ri < V.nr ? { ri: ri } : null;")

# control bar of the window: 'Rows by tree' puts the rows in the order of the tree (only shown while the rows are in another order)
rep('<span class="lbl" id="tdColLbl" hidden>Tips</span>',
    '<button class="chip" id="tdRowsBy" hidden style="font-size:11px;padding:0 8px" title="Put the rows in the order of this tree (the lines become horizontal)">Rows by tree</button><span class="lbl" id="tdColLbl" hidden>Tips</span>')
rep('$("tdClose").addEventListener("click", function () { tdToggle(false); });', '$("tdClose").addEventListener("click", function () { tdToggle(false); }); $("tdRowsBy").addEventListener("click", function () { tdSyncOrder(); });')

# opening the window keeps a name (or curation) order
rep('if (TD.on && S.order !== "ml:" + TD.set) tdSyncOrder(); else { sizeScroller(); schedule(); }', 'if (TD.on && S.order !== "ml:" + TD.set && !tdKeepOrder()) tdSyncOrder(); else { sizeScroller(); schedule(); }')
rep('if (TD.on && S.order !== "ml:" + TD.set) { S.order = "ml:" + TD.set; $("orderSel").value = S.order; }', 'if (TD.on && S.order !== "ml:" + TD.set && !tdKeepOrder()) { S.order = "ml:" + TD.set; $("orderSel").value = S.order; }')
rep("(rows take the order of the tree)\">Tree window", "(rows take the order of the tree; with Order = Name the tree keeps its own order and a line joins each tip to its row)\">Tree window")

open(P, "w", encoding="utf-8", newline="").write(t)
print("patched", P, len(t), "bytes")
