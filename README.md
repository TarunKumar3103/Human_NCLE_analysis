# Human NCLE–NED association analysis

Tests whether human proteins with a native non-covalent lasso entanglement (NCLE) are more likely to show non-exponential degradation (NED vs ED in McShane et al., *Cell* 2016, RPE-1 cells). The repo contains the recovered original analysis code, frozen inputs, the reviewed reference workbook, reproduction and audit scripts, and an independent verification.

## Results

| Analysis | Model | N | Term | OR | 95% CI | P |
|---|---|---|---|---|---|---|
| Primary: NED vs ED | Crude | 705 | NCLE | 0.96 | 0.64–1.44 | 0.835 |
| | Model 1 | 705 | NCLE | 0.95 | 0.63–1.44 | 0.818 |
| | Model 2 | 528 | NCLE | 0.93 | 0.60–1.44 | 0.741 |
| | Model 2 | 528 | Complex member | 9.16 | 1.23–68.24 | 0.031 |
| Sensitivity: NED vs ED + UN | Crude | 1,191 | NCLE | 0.93 | 0.63–1.38 | 0.763 |
| | Model 1 | 1,191 | NCLE | 0.95 | 0.64–1.41 | 0.800 |
| | Model 2 | 901 | NCLE | 0.90 | 0.60–1.36 | 0.623 |
| | Model 2 | 901 | Complex member | 6.73 | 0.91–49.58 | 0.061 |

No NCLE–NED association was detected in either cohort or any model; every NCLE estimate is close to 1, and every NCLE CI includes 1. This is not proof of no effect. The complex-membership estimate is compatible with the complex-membership hypothesis but imprecise; see Limitations.

**Methods in brief.**

- **Crude:** sample odds ratio, Woolf log-OR 95% CI, two-sided Fisher exact P.
- **Model 1:** logistic regression `NED ~ NCLE + z(protein length)`.
- **Model 2:** adds `Complex_member` (1 = complex/oligomer member, 0 = explicit monomer, unknown = missing and excluded).
- **Length:** standardized once within each full cohort (sample SD).
- **Regression CIs and P:** Wald.

Cohort filters, in order:

- average pLDDT ≥ 70
- binary NCLE (any covalent lasso record excluded)
- not knotted
- FASTA sequence available
- not a membrane protein
- not a secreted protein

Full details are in `README.txt` and `documentation/`.

## Reproduce

Python 3.12 and the pinned packages in `requirements.txt` (statsmodels is needed only for the independent check). Everything runs offline from local files.

~~~
python -m pip install -r requirements.txt statsmodels

python scripts/export_complex_annotations.py
python scripts/reproduce_results.py
python scripts/audit_raw_sources.py
python scripts/audit_disorder.py
python scripts/prepare_mobidb_ids.py
python scripts/independent_checks.py --human .
~~~

Expected output:

- `export_complex_annotations.py`: `Recorded-decision reconstruction: PASS`. That is 1,191 evidence records and 205 reviewed, giving 28 provisional complex calls and 177 unresolved.
- `reproduce_results.py`: `Numerical reproduction: PASS`. All 12 reported rows match, and the six Model 2 slope rows agree with a second optimizer.
- `audit_raw_sources.py`: every cohort-flow row is `True`. That is 1043 → 1034 → 1023 → 1022 → 725 for primary and 1836 → 1810 → 1793 → 1792 → 1230 for sensitivity, through membrane exclusion.
- `audit_disorder.py`: 13,690 MobiDB rows, 12,609 fractions, 18 coordinate warnings.
- `prepare_mobidb_ids.py`: 13,690 unique IDs.
- `independent_checks.py`: **19/19 independent checks passed**.

What `independent_checks.py` checks, without importing any package code:

- It rebuilds both cohorts from the raw files.
- It confirms that the final 705 and 1,191 proteins are subsets of the rebuilt pre-secreted cohorts, with matching NCLE, length, pLDDT and degradation profile.
- It refits all 12 result rows with statsmodels.

## Layout

| Folder | Contents |
|---|---|
| `scripts/` | `statistical_analysis.py` (recovered original, unchanged), reproduction runner, audit utilities, MobiDB query script, independent check |
| `data/` | Reviewed reference workbook, processed MobiDB export, exported complex annotations and evidence snapshot |
| `raw_sources/` | McShane human supplement, Jiang NCLE and pLDDT tables, AF4 knots, UniProt proteome FASTA, membrane-keyword FASTA |
| `results/` | Validated outputs (frozen replay, raw-source audit, disorder audit, MobiDB query inputs) |
| `historical_audits/` | Older audits, repository notebook, Quyen's original disorder script; their cutoffs are not the current protocol |
| `documentation/` | Reproduction guide (PDF/DOCX), code compendium, independent-review prompt |

`SOURCE_MANIFEST.json` lists file hashes and provenance.

## Limitations

- **Secreted filter not independently replayed:** the historical secreted-keyword export was not recovered. Every step up to and including membrane exclusion reproduces from raw sources. The final cohorts are exact subsets of the rebuilt pre-secreted cohorts: 20 proteins removed in primary, 39 in sensitivity.
- **Complex annotations:** the original text classifier was not recovered. Saved UniProt evidence and manual-review decisions reproduce the labels exactly, but 28 calls are provisional and 177 proteins stay unresolved (missing).
- **Complex-membership term:** only one NED protein appears among the explicit monomers. The Wald P is conservative here; likelihood-ratio P is 0.002 (primary) and 0.010 (sensitivity). Annotation availability is itself associated with NED, so selection bias is possible.
- **Isoform IDs:** exact accession matching drops 353 McShane ED/NED isoform IDs. Mapping them to canonical accessions would add 148 proteins before the secreted filter. The crude OR then becomes 0.87 (P = 0.50), so the null result holds.
- **MobiDB disorder** (not used in any statistics):
  - 33 fractions use a current UniProt length that differs from the proteome FASTA the coordinates match; Q8NI36 is one.
  - Proteins with no selected annotation are probably close to 0% disorder rather than missing.
- **Interpretation:** the ORs are associations, not causal effects.

## Sources

McShane et al., *Cell* 167, 803–815 (2016), https://doi.org/10.1016/j.cell.2016.09.015 · MobiDB, https://mobidb.org · UniProt, https://www.uniprot.org
