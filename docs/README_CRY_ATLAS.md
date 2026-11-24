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

## Trees, ancestors, results of 2026-11-24 (added after the first publish)
* **Trees** (`scripts/20_trees.py`): IQ-TREE 2.3.6, fixed LG+G4, `--fast` search, SH-aLRT support (1000), seed 20261124. Full-toxin trees use the 1,402 (structure frame) / 1,261 (sequence frame) columns occupied by >= 20 % of the entries; domain trees (DI, DII, DIII) use the columns of the module and the entries with >= 100 (DI, DII) or 70 (DIII) residues of it. Eight trees in `site/data/cry/trees.json`; the tanglegram opens on the full toxin against Domain III.
* **Ancestors** (`scripts/asr/`, copy of the NODE-MSA Atlas ASR tools; PAML 4.10.10 codeml LG+G4 with the IQ-TREE branch lengths + GRASP 2.0.0): 669 internal nodes for each of the two full-toxin trees (`asr_struct_A.json`, `asr_seq_A.json`, index `asr_index.json`). Median node: mean posterior probability 0.999, PAML = GRASP at 100 % of the residues, 99.5 % identical to the nearest extant toxin; the root is uncertain (mean PP 0.60).
  GRASP's default indel method (BEP) did not finish in the 2-hour limit (killed after 2 h), PSP with all 1,402 columns ran out of Java heap at 6 GB, and PSP with CONTIGUOUS column blocks failed ("Failed to find optimal path for ancestor") because whole clades (Cry2, Cry3 ...) lack the protoxin half. The final run uses PSP with **interleaved** blocks (columns c, c + 4, c + 8 ...): `asr_parse.run_grasp` was changed accordingly. Environment: `ASR_MIN_OCC=0.2` (the occupancy of the tree columns), `ASR_GRASP_CHUNKS=4`, `ASR_GRASP_MEM=6g`, `py asr_export.py cry <key> --indel PSP`. PAML (about 3.5 h with the CPU shared) is cached in `scripts/asr/out/` (not in the repository).
* **Not done**: species-tree window (almost all hosts are one species, Bacillus thuringiensis); no ancestors for the domain trees (the ASR tools cut the alignment by subunit, not by module).
