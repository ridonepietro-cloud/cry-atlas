# Cry Atlas (insecticidal 3-domain proteins of the BPPRC database)

Separate Atlas artifact built from the NODE-MSA Atlas pipeline (page copy + pack builders) for the **Cry class of bpprc-db.org**
(`https://www.bpprc-db.org/categorize_database_cry`, downloaded 2026-11-24). Nothing of the NODE-MSA Atlas (`atlas/site`, its packs, its published artifact) was edited.

## Data sources (all public)
| what | source | file |
|---|---|---|
| Cry names, old / other names, accession, year (808) | BPPRC category page `categorize_database_cry` | `data/bpprc_cry_table.tsv` |
| full-length sequences (797 records, 795 names) | BPPRC `category_download` (category = cry) | `raw/fasta/cry_full.fasta` |
| BPPRC domain sequences D1 / D2 / D3, holotype list | BPPRC `threedomain_download`, `holotype_download` | `raw/fasta/` |
| structures table (PDB ids per BPPRC name) | BPPRC `structures/` | `data/bpprc_cry_structures.tsv` |
| specificity (target order / species, active yes / no) | BPPRC `export_association` (curator column dropped) | `work/specificity.json` |
| organism, strain, plasmid, coded_by, NCBI sequence | NCBI GenPept (efetch) | `work/ncbi/genpept.json` |
| UniProt entries of identical sequences, AlphaFold DB models | UniParc (MD5) -> UniProtKB -> AlphaFold DB API; model accepted only when its sequence equals the row | `work/afdb/map.json` |
| PDB entries | RCSB (accession search + full-text candidates, entity sequences by GraphQL, local matching >= 90 %) | `work/pdb/matches.json` |
| genomic context | NCBI IPG + nuccore GenBank windows (gene +- 15 kb) | `site/data/cry/gc.json`, `work/gc_windows/` |

## Pipeline (scripts in `scripts/`, run with `py -I`; run order)
01 parse BPPRC table -> 02 FASTA + domains -> 03 clean / merge exact duplicates -> 04 GenPept -> 05 compare BPPRC vs NCBI sequences -> 06 BPPRC structures table ->
07 UniProt + AlphaFold DB map -> 08b/08c PDB candidates and entity matching -> 09 rows -> 10 / 10b sequence frame (first pass, then segment-wise MAFFT L-INS-i) -> 15 modules ->
11 specificity -> 12 download structures -> 13 RING 4.0 -> 14 assemble rows (models) -> 16b structure frame (segment-wise FoldMason) -> 17 classes -> 19 build the pack (`cry_pack.py`) ->
`verify_cry_pack.py` -> 18 / 18a / 18b genome map and synteny -> 20 trees (IQ-TREE) -> `asr/` ancestors (PAML + GRASP) -> `patch_cry_page.py` (page) -> publish.
Large raw files (AlphaFold DB, PDB, RING edge files) are in `%TEMP%\cry_raw` (re-downloadable; manifest `work/structures_manifest.tsv` with URL, size and SHA-256).

## Decisions (and why)
* **Rows** = distinct BPPRC sequences (670 from 795 names); exact duplicates are merged and all their names listed in the row. The BPPRC sequence is the entity (NCBI differences are listed in the row info).
* **Primary model** of a row = the model that places most residues of its sequence (full-length first; experimental wins a tie); the others are alternative models of the row.
* **Modules** (N-term | DI | DII | DIII | C-term, the C-terminal protoxin half) = consensus of the BPPRC domain annotation over 574 annotated rows, carried by alignment so every row gets the same boundaries.
* **NODE** of a monomeric toxin = the packages' rule applied to contacts that join two modules (domains) inside the chain: a column is strict NODE when its inter-module contacts over all models have >= 2 interaction types (one not VDW) and occur in >= 5 models.
* **Frames**: both alignments are built per module and joined (a single-pass alignment of full-length toxins stretched DI to 610 columns): structure-based = FoldMason of the AlphaFold DB models of each module (+ MAFFT --add for rows without an exact model); sequence-based = MAFFT L-INS-i of each module.
