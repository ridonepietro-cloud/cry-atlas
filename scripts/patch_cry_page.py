"""Page of the Cry Atlas: patch of a COPY of the NODE-MSA Atlas page (site/base/index_atlas_20261122.html) -> site/index.html.
One system / one pack ('cry'): branding, tab bar, pack table, fallbacks, CSS that hides the FIT-specific windows, wording for a single-chain toxin with domains.
Idempotent (always starts from the base copy).  Run: py -I patch_cry_page.py   then (separate step) the download-buttons patch of the Atlas page:  py -I <atlas>/patch_atlas_downloads_20261124.py site/index.html"""
import os, sys, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cry_common import SITE
BASE = os.path.join(SITE, "base", "index_atlas_20261122.html"); OUT = os.path.join(SITE, "index.html")
t = open(BASE, encoding="utf-8", newline="").read()
t = t.replace(chr(13) + chr(10), chr(10))                 # the base page has CRLF line ends; the Cry page is written with LF


def rep(old, new, count=1):
    global t
    assert t.count(old) == count, (t.count(old), old[:90])
    t = t.replace(old, new)


def rep_re(pat, new, flags=re.S):
    global t
    n = len(re.findall(pat, t, flags)); assert n == 1, (n, pat[:80])
    t = re.sub(pat, lambda m: new, t, flags=flags)


rep("<title>NODE-MSA Atlas</title>", "<title>Cry Atlas</title>")
rep("<h1>NODE-MSA Atlas</h1>", "<h1>Cry Atlas <span class=\"stamp\" style=\"font-weight:400\">3-domain insecticidal proteins of the BPPRC database</span></h1>")
rep_re(r'<nav class="tabs" role="tablist" aria-label="System" id="sysTabs">.*?</nav>',
       '<nav class="tabs" role="tablist" aria-label="System" id="sysTabs">\n      <button role="tab" aria-selected="true" data-sys="cry" title="Cry toxins: every sequence of the Cry class of bpprc-db.org with its structure models, domains, contacts, trees and genome maps">Cry toxins</button>\n    </nav>')
rep_re(r'var PACKS = \{.*?\};\nvar SYSNAMES = \[[^\]]*\], DEFPACK = \{[^}]*\};',
       'var PACKS = { cry: { sys: "cry", dir: "cry", label: "Cry toxins (BPPRC Cry class)" } };\nvar SYSNAMES = ["cry"], DEFPACK = {};')
rep('if (SYSNAMES.indexOf(h) >= 0) return packOf(h) || "fit";', 'if (SYSNAMES.indexOf(h) >= 0) return packOf(h) || "cry";')
rep('return (SYSNAMES.indexOf(s) >= 0 && packOf(s)) || "fit";', 'return (SYSNAMES.indexOf(s) >= 0 && packOf(s)) || "cry";')
rep('''    document.querySelectorAll('[data-scope="inter"]').forEach(function (b) { b.textContent = "Between protomers"; });
''', '''    document.querySelectorAll('[data-scope="inter"]').forEach(function (b) { b.textContent = "Between protomers"; });
    if (SYS === "cry") { document.querySelectorAll('[data-scope="intra"]').forEach(function (b) { b.textContent = "Contacts"; b.title = "All stored contacts of the toxin chain (side-chain contacts, and every contact that joins two domains)"; }); }
''')
# ---- 'Colour entries by' bar (the features framework of the IRM tab) also for the Cry tab: attributes of the rows from data/cry/features.json
rep('function ftOn() { return SYS === "irm" && !!(S.meta && S.meta.irm && S.meta.irm.features); }',
    'function ftOn() { var ft0 = S.meta && (S.meta.irm || S.meta.cry); return (SYS === "irm" || SYS === "cry") && !!(ft0 && ft0.features); }')
rep('fetchJSON("data/" + PACK + "/" + S.meta.irm.features.file, 1)', 'fetchJSON("data/" + PACK + "/" + (S.meta.irm || S.meta.cry).features.file, 1)')
rep('if (!FT.data || SYS !== "irm") return false;', 'if (!FT.data || (SYS !== "irm" && SYS !== "cry")) return false;')
rep('''  h += '<button class="chip" id="ftOpen" aria-pressed="' + FT.open + '" title="Core architectures, group × architecture matrix, appendage classes and the association between them, with the caveats">Architectures</button>';''',
    '''  if (SYS === "irm") h += '<button class="chip" id="ftOpen" aria-pressed="' + FT.open + '" title="Core architectures, group × architecture matrix, appendage classes and the association between them, with the caveats">Architectures</button>';''')
rep('''  $("ftOpen").addEventListener("click", function () { ftSetOpen(!FT.open); });
}''', '''  if ($("ftOpen")) $("ftOpen").addEventListener("click", function () { ftSetOpen(!FT.open); });
}''')
rep('var h = "<section><h2>Core architecture and appendages</h2>', 'var h = "<section><h2>" + (SYS === "cry" ? "Attributes of the entry" : "Core architecture and appendages") + "</h2>')
rep('function lkEntryHTML(r) {', 'function lkEntryHTML(r) {\n  if (SYS === "cry") return ftEntryHTML(r);')
# ---- tanglegram: the right-hand tree follows the alignment frame on screen (struct_DIII <-> seq_DIII); the colour option is not about stators
rep('t.left = S.trees.trees[f + "_A"] ? f + "_A" : td0[0]; var alt = td0[1];', 't.left = S.trees.trees[f + "_A"] ? f + "_A" : td0[0]; var alt = td0[1] && (S.trees.trees[f + td0[1].slice(td0[1].indexOf("_"))] ? f + td0[1].slice(td0[1].indexOf("_")) : td0[1]);')
rep('>by stator class</option>', '>by family class</option>')
css = ("html[data-sys=\"cry\"] #fpToggle, html[data-sys=\"cry\"] #fpPanel, html[data-sys=\"cry\"] [data-c=\"hyd\"], html[data-sys=\"cry\"] [data-c=\"hydw\"], html[data-sys=\"cry\"] [data-c=\"hydx\"], "
       "html[data-sys=\"cry\"] #stBtn, html[data-sys=\"cry\"] [data-scope=\"inter\"], html[data-sys=\"cry\"] [data-layer=\"plug\"], html[data-sys=\"cry\"] #plugFilter { display: none !important; }\n")
rep("</style>\n\n<div id=\"app\">", css + "</style>\n\n<div id=\"app\">")
open(OUT, "w", encoding="utf-8", newline="").write(t)
print("written", OUT, len(t), "bytes")
