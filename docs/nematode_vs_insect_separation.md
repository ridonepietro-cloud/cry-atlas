# Nematode-active versus insect-active Cry toxins: which columns and contacts separate them (2026-11-27)

Groups from the BPPRC specificity export (per sequence, union over its names): {'Insect-active': 173, 'No data': 461, 'Other organism only': 1, 'Nematode-active': 19, 'Tested, none active': 16}.
Positives = 17 nematode-active toxins of six families (Cry5, Cry12, Cry13, Cry14, Cry21, Cry95); negatives = 173 insect-active toxins of 28 families.

**Alignment columns** (structure-based frame, columns with >= 30 % occupancy): 7 residue-column features pass all criteria (>= 80 % of the positives, <= 10 % of the negatives, in the majority of the rows of >= 5 of the 6 nematode families).
Leave-one-nematode-family-out AUC per held-out family: {'Cry95': 0.41, 'Cry5': 1.0, 'Cry21': 1.0, 'Cry13': 0.95, 'Cry14': 1.0, 'Cry12': 1.0}, mean 0.89; random pseudo-nematode families (150 draws): mean 0.52, 95th percentile 0.75, p = 0.013.

**Contacts** (stored contacts of the primary models, column pairs present in >= 3 % of the rows): 10 column pairs pass the same criteria. Leave-one-family-out mean AUC 0.95; random families mean 0.52, 95th percentile 0.71, p = 0.007.
Most of the passing contacts lie in the C-terminal (protoxin) module: a hint that the nematode toxins share a protoxin architecture that the insect-active toxins in the sample do not share at the same aligned positions.

Caveats: 17 positives in 6 families (Cry12, Cry13 and Cry95 are single entries); the nematode toxins are close relatives, so a signature is first a lineage marker; the family-label permutation test controls for this only in part; the activity labels are BPPRC records, which are incomplete (461 of 670 sequences have no specificity record).
Files: nematode_signature_columns.tsv, nematode_signature_contacts.tsv; scripts 24_groups_and_separation.py, 27_apply_layers.py.
