# Cry Atlas: quick guide to selecting, navigating and visualising entries

Open **https://ridonepietro-cloud.github.io/cry-atlas/site/** in a desktop browser (Chrome, Edge or Firefox). The first load takes a few seconds (about 40 MB of data). The page holds the 670 distinct BPPRC Cry sequences (795 names; entries with identical sequences share one row).

## 1. Pick the entries you care about

Everything below acts on the entries currently *visible*. Combine the filters freely.

| What you want | How |
|---|---|
| One family or a few | Click the **family chips** (Cry1 · Cry3/7/8 · Cry2/11/18… · Cry21/Cry5/Cry14… · …). Click again to hide or show a class. |
| One kind of entry (nematode-active, domain swap, Lepidoptera-active, one year…) | **Colour entries by** → pick an attribute → click the **value chips** (e.g. "Insect-active 173") to hide or show those entries. The number on each chip is the entry count. Attributes: Activity group, Domain trees (does one domain differ?), Target order, Source organism, Year, Length, Structure model. "reset filters" clears them. |
| One name | Type in **Filter entries** (e.g. `Cry5B`, `Cry1Ac`). |
| A hand-picked set | Tick the boxes at the left of the entry names. |
| Everything visible, or only entries with a structure | **Selection → Select…** → *All visible*, *With an experimental complex*, *With any structure*, *Invert visible*, *Clear*. |

**Open as subset** (needs at least 2 ticked entries) puts them in their own workspace chip beside **All entries (670)**. Columns that are empty in those entries are hidden; it is a projection of the full alignment, not a new alignment. Click **All entries (670)** to go back. **Export list** saves the entry list.

## 2. Look at the alignment

* **View → Alignment**; **Structure** (FoldMason on AlphaFold models, per module) or **Sequence** (MAFFT L-INS-i, per module). The modules are N-term | Domain I | Domain II | Domain III | C-term (protoxin half).
* **Order**: *Tree order* (guide tree of the Domain II alignment, the default), *ML tree* (IQ-TREE of the frame on screen) or *Name*.
* **Colour**: Confidence (pLDDT), Module, Occupancy, Residue type. **Zoom / Rows** sliders, **Fit**, **Fit all**.
* **Overlays**: *NODE* (columns where the domains touch each other), *Contact dots*, *Galactose sites* (the marked sites, see §5), *Column strip*; **Choose layers** lists them all. The **Partners** chips (H-bond, ionic, π-stack…) filter the contact types.
* **↓ FASTA** downloads the sequences of the visible entries (the arrow next to it opens the format options).

## 3. Open one entry

* **Click an entry name**: its structure appears (Modules, Confidence, N→C colouring, Reset view, Copy FASTA) with the **entry card**: all BPPRC names, accession, organism, year, specificity records, domain boundaries, activity group, domain-tree call.
* **Click a residue** in the alignment: its contact partners are listed.

## 4. Trees

* **View → Trees**. Choose the tree on screen: full toxin or Domain I / II / III, each on the structure or the sequence frame. Click a node to pick that clade (alt-click removes). **↓ Newick** exports the tree on screen, **Pop out** detaches it.
* **Tree window** (next to Alignment) draws the tree beside the alignment rows.
* **Tanglegram of two trees**: Left / Right trees of your choice (default Domain II against Domain III). **by family class** colours the lines by family; **by discordance** highlights the entries whose position differs between the two trees, which is the domain-swap view.

## 5. Other windows

* **Ancestors**: reconstructed ancestors of every node of the 8 trees (Consensus, GRASP, PAML), with posterior probability per site.
* **Genome map**: tick entries first; shows ±5/10/15 kb around each gene with synteny between the ticked entries.
* **Species tree** (button at the top, next to *Cry toxins*): GTDB tree of the source species with one ring per Cry class; click a box to list its entries.
* **Marked sites** (diamonds and stars): diamonds = residues within 4.5 Å of the N-acetylgalactosamine in 4ARY (Cry1Ac); stars = the seven alignment columns that separate the nematode families from insect-active toxins. Tick *Galactose sites* in Overlays; this one toggle switches both the diamonds and the stars.

## 6. Three ready-made routes

1. **Nematode-active vs insect-active.** Colour entries by → *Activity group* → click *Insect-active*, *Other organism only*, *Tested, none active* and *No data* to hide them (19 nematode-active entries remain) → Select… → *All visible* → **Open as subset**. Show *Insect-active* again to compare. Switch *Overlays* to the marked sites to see the stars.
2. **Domain swaps.** Colour entries by → *Domain trees* → hide *concordant*, *not assessable* and *not in all domain trees*: 11 entries remain whose neighbours in one domain tree disagree with the other two (8 are Domain III in Cry1B/C/E/J; 2 differ in all three). Click a name for the entry card, or open the **Tanglegram** with *by discordance*.
3. **The galactose site.** Overlays → *Galactose sites*, Structure frame, Colour: Residue type; zoom to Domain III. The eight diamond columns are the 4ARY contact residues, so what every entry carries at those columns is visible row by row (conservation is low: 7-27 % of entries share the Cry1Ac residue).

## Limits to keep in mind

* This public copy is static: *Save realign job*, comments and *Alignment curation* only work in the private Claude artifact version. The filters, subsets, trees, ancestors and downloads work here.
* The nematode signature and the domain-tree calls are **candidates**: 17 positives in 6 families (the signature may be a lineage marker), and 320 of 667 entries cannot be judged for domain swaps because their nearest relatives do not share a name group.
* 70 entries have no structure model; 49 use a flagged near-identical donor model. Activity data exist for only 209 of 670 sequences, and 10 of the 20 nematode-active names were tested on one species.
* Notes and tables behind each analysis are in the repository's `docs/` folder.
