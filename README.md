# Cry Atlas

An interactive atlas of the **670 distinct Cry toxin sequences** (795 names) of the Cry class of the
[Bacterial Pesticidal Protein Resource Center](https://www.bpprc-db.org/) (BPPRC), downloaded 2026-11-24.

**Live page:** https://ridonepietro-cloud.github.io/cry-atlas/ (opens `site/index.html`)

What it shows, per toxin: the sequence in a structure-based and a sequence-based alignment (both built domain by domain: N-terminus, Domain I, II, III,
C-terminal protoxin half), AlphaFold DB and PDB models, RING 4.0 contacts and the residues that connect domains (NODE residues), maximum-likelihood trees
of the whole toxin and of each domain (with a tanglegram to compare them), genome maps with gene synteny, BPPRC target-order data, and the BPPRC names, years and holotype flags.

The page is a derivative of the NODE-MSA Atlas viewer (3Dmol.js for the structures). It was made for the claude.ai artifact runtime; as a static site on GitHub Pages the
functions that need that runtime (saved subsets, curation, comments, the save-file dialog) are inactive. Everything else runs in the browser.

## Layout
| folder | content |
|---|---|
| `site/` | the page (`index.html`) and its data pack `site/data/cry/` (rows, alignments, contacts, structure bundles, trees, genome maps, ancestors when present) |
| `scripts/` | the numbered pipeline (see `docs/README_CRY_ATLAS.md` for the run order and the decisions), the pack writer `cry_pack.py`, the verifier `verify_cry_pack.py`, `asr/` (PAML + GRASP wrappers) |
| `data/` | tables parsed from the BPPRC pages (names, accessions, domains, structures) |
| `docs/` | pipeline notes |

Raw downloads (AlphaFold DB and PDB files, RING edge files, PAML / GRASP run folders) are not in the repository; the scripts re-create them. The page that
`scripts/patch_cry_page.py` patches (`site/base/index_atlas_20261122.html`) is a copy of the NODE-MSA Atlas page and is not included.

## Data sources and terms
* **BPPRC** (bpprc-db.org): names, sequences, domain boundaries, holotypes, structure table, specificity data (manually curated from papers and patents, each source cited there).
  Please cite BPPRC as described at https://www.bpprc-db.org/cite/ and check the original publication of any specificity record you rely on.
* **AlphaFold DB** models (CC BY 4.0), **RCSB PDB** entries (CC0), **UniProt / UniParc** (CC BY 4.0), **NCBI** GenPept and nucleotide records (public data, NCBI usage policy).
* Tools: MAFFT, FoldMason, IQ-TREE, PAML (codeml), GRASP, RING 4.0 (academic licence; only the contact tables derived from it are in this repository), 3Dmol.js.
* Alignments, trees, contacts and annotations here are computed from those sources and may contain errors; sequence differences between BPPRC and NCBI are listed in the row notes.

Code in `scripts/` and `site/index.html` has no licence file yet; ask before reusing it.
