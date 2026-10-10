# Independent reproducibility audit: human NCLE/decay data curation

Audit date: 2026-09-26  
Scope: data curation and cohort construction only. No association testing, matching, simulation, or secretory-protein exclusion was performed.

## 1. Verdict

**PASS WITH WARNINGS.**

The central cohort logic is reproducible from the supplied materials:

- The pLDDT-qualified master contains **13,690 unique proteins**.
- Exact-accession integration produces **1,043 ED/NED proteins**.
- Applying the documented covalent-lasso, knot, FASTA, and membrane filters produces **725 non-membrane core proteins**.
- Requiring a nonmissing MobiDB disorder fraction produces **606 proteins**.
- The four ED/NED × entanglement group counts exactly reproduce the reported values.
- Six fields that can be independently reconstructed from McShane and Jiang sources match the processed 13,690-row export for every protein: mean pLDDT, pLDDT SD, entanglement value, entanglement label, ED/NED label, and binary decay outcome.

The 606-protein membership is therefore reproducible. It is **conditionally ready for simulation candidate selection**, but not yet a fully archived/final deliverable, because the three generated Excel workbooks were not supplied, the production MobiDB script and raw/full query result were not supplied, the MobiDB method description conflicts with the attached code, and one protein in the 606 cohort has incompatible displayed disorder coordinates.

## 2. Files received and missing

### Received and usable

- McShane human degradation workbook, sheet `FinalTable_RPE_`.
- Jiang mean pLDDT file: `Human_Avg_pLDDTs.csv`.
- Jiang entanglement file: `Human_AF_combined_20250614.csv`.
- AlphaKnot file: `AF4_knots.dat`.
- Human reference proteome FASTA: `UP000005640_9606.fasta`.
- UniProt membrane-keyword FASTA.
- A 13,690-row processed broad-master CSV export embedded in the supplied text file. This contains MobiDB result columns.
- Quyen's original MobiDB README, query script, and five-protein example output.

### Missing or incomplete

- `Human_master_pLDDT_ge70_with_MobiDB_and_membrane_QC.xlsx`.
- `Human_ED_NED_nonmembrane_core_pLDDT_ge70.xlsx`.
- `Human_ED_NED_nonmembrane_MobiDB_complete.xlsx`.
- The full production MobiDB TSV used to construct the master.
- The production query script that generated `query_status` and `annotation_tracks`. Quyen's supplied script does not emit those fields.
- Raw MobiDB JSON responses or another immutable API snapshot.
- Both supplied repository ZIP copies lack a valid ZIP central directory. The relevant payloads could be recovered and parsed, but their integrity cannot be checked against the original GitHub archive without an external reference hash.

Exact paths, byte sizes, and SHA-256 hashes are in `results/file_inventory.csv`.

## 3. Schemas and identifier rules used

| Source | Identifier/field | Audit rule |
| --- | --- | --- |
| McShane | `Protein IDs (Uniprot)` | Trim whitespace; otherwise use exact accession. No silent canonicalization. |
| McShane | `Degradation profile` | `ED → 0`, `NED → 1`; `UN` is not included in the primary ED/NED cohort. |
| Jiang pLDDT | `gene`, `<pLDDT>`, `std` | Exact `gene` accession; retain mean pLDDT **≥ 70**. No protein has exactly 70. |
| Jiang entanglement | `gene`, `ENT-ID`, `CCBond` | Group by accession. Any `CCBond=True` excludes the protein as a covalent lasso. Otherwise at least one `CCBond=False` means native NCLE; no entanglement rows means no detected NCLE. |
| AlphaKnot | `Uniprot`, `Organism` | Require organism exactly `Human`; remove AlphaFold fragment suffix `-F#`; exact accession match thereafter. |
| Human FASTA | second pipe-delimited token | Exact accession match. |
| Membrane FASTA | second pipe-delimited token | Exact accession match; keyword `Membrane` is used as the approved operational definition. |
| MobiDB | processed master columns | Complete when `Fraction_IDP_IDR` is nonmissing; in the supplied export this is exactly equivalent to `Disorder_query_status == "ok"`. |

The entanglement source supports the aggregation rule operationally: all **9,247** proteins classified as not entangled have no `ENT-ID` and no `CCBond` value; **10,768** have one or more noncovalent records; **573** have at least one covalent record and are excluded. The label should nevertheless be interpreted as **“no detected NCLE in this dataset,” not proof that entanglement is biologically impossible**.

## 4. Cohort flow

| Sequential step | Proteins remaining | Removed at this step |
| --- | ---: | ---: |
| Mean pLDDT ≥70 broad master | 13,690 | — |
| Retain ED/NED | 1,043 | 12,647 |
| Exclude covalent-lasso entanglements | 1,034 | 9 |
| Exclude knotted proteins | 1,023 | 11 |
| Require human proteome FASTA accession | 1,022 | 1 |
| Exclude UniProt membrane-keyword proteins | 725 | 297 |
| Require nonmissing MobiDB fraction | **606** | 119 |

Secretory proteins were not removed in this reconstruction. No IDP/IDR-fraction cutoff was applied.

## 5. Group counts and MobiDB availability

| Analysis group | Core | MobiDB complete | Missing MobiDB | Missing % |
| --- | ---: | ---: | ---: | ---: |
| ED, entangled | 415 | 345 | 70 | 16.87% |
| ED, not entangled | 169 | 138 | 31 | 18.34% |
| NED, entangled | 98 | 84 | 14 | 14.29% |
| NED, not entangled | 43 | 39 | 4 | 9.30% |
| **Total** | **725** | **606** | **119** | **16.41%** |

These are descriptive QC counts only. No odds ratios, P-values, or regression models were calculated.

## 6. Benchmark comparison

All expected/reported counts were reproduced exactly except one:

- McShane: 3,133 total; 1,474 ED; 296 NED; 1,363 UN; 1,770 ED/NED.
- Jiang pLDDT: 20,588 unique proteins; 13,690 with mean pLDDT ≥70.
- Broad entanglement: 8,416 native-entangled; 4,782 not-entangled; 492 covalent-lasso excluded.
- MobiDB status: 12,609 `ok`; 1,042 no selected annotation; 31 no selected annotation plus missing length; 8 missing length.
- FASTA availability: 13,594 of 13,690.
- Cohorts: 1,043 ED/NED; 1,022 FASTA-complete pre-membrane; 725 non-membrane; 606 MobiDB-complete.
- **Knot discrepancy:** the earlier notebook transcript reported 123 knotted proteins in the broad master, but the supplied/recovered `AF4_knots.dat` produces **122** exact matches. The downstream ED/NED cohort contains 11 knot exclusions and still reproduces 1,023 and 1,022 exactly. This is therefore a broad-master provenance/count discrepancy, not a cause of the 606-cohort difference.

See `results/benchmark_comparison.csv` for every benchmark.

## 7. Generated-workbook comparison

The exact three Excel files named above were not present, so their sheet schemas, row ordering, formulas, cell values, and ID sets could not be compared directly.

The supplied 13,690-row processed export was compared with a fresh reconstruction from McShane/Jiang sources:

| Field comparison | Mismatches |
| --- | ---: |
| ID set | 0 |
| Mean pLDDT | 0 |
| pLDDT SD | 0 |
| Entanglement binary value | 0 |
| Entanglement text status | 0 |
| ED/NED profile | 0 |
| Binary decay outcome | 0 |

This validates the broad integration logic, but it is not a substitute for workbook-level validation. The MobiDB columns can only be checked internally because their production TSV/script/raw responses were not supplied.

## 8. Duplicate, missingness, range, and logic checks

- No duplicate or missing IDs occur in McShane, pLDDT, the processed master, either FASTA, or human AlphaKnot records after the documented fragment normalization.
- All 13,690 master IDs are unique and the `Index` is sequential 1–13,690.
- No retained mean pLDDT is below 70. The minimum retained value is 70.000377; the maximum excluded value is 69.997463.
- All nonmissing disorder fractions fall in [0,1].
- Every fraction equals `Disorder_residues / MobiDB_length` to the saved eight-decimal precision.
- No stored disorder-residue count exceeds the stored MobiDB length.
- Every `ok` record has a disorder fraction; every non-`ok` record has a missing fraction.
- All 119 missing disorder values in the 725-protein core are `mobidb_no_selected_annotation`, not numerical zeros.
- The 13,690-row broad master intentionally blanks `UN` along with proteins absent from McShane. There are 793 pLDDT-qualified UN proteins in the source. Therefore the broad master alone cannot support the planned `NED` versus `ED + UN` sensitivity analysis without re-merging the original McShane profile.

## 9. MobiDB findings

### 9.1 Fraction arithmetic passes

The saved fraction calculation is correct to eight decimal places for all 12,609 complete records.

### 9.2 Eighteen displayed interval/count inconsistencies

Eighteen broad-master proteins have disorder intervals that extend past the saved UniProt length. Summing the displayed intervals gives a larger number than `Disorder_residues`, but clipping the intervals to positions 1 through the saved length reproduces the stored count in all 18 cases. This is consistent with sequence-coordinate/version mismatch, not a random arithmetic error.

Only **Q8NI36** is in the final 606:

- saved length: 895;
- displayed intervals: `1-63;722-734;807-809;949-951`;
- un-clipped displayed count: 82 residues;
- count clipped to 1–895: 79 residues;
- saved fraction: 79/895 = 0.08826816.

Recommendation: preserve both raw and clipped intervals in separate fields, add an explicit coordinate-version/QC flag, and requery Q8NI36 against the same sequence version used for the structure. Do not silently change the fraction or delete the protein.

### 9.3 README/code disagreement

The README says the script selects the highest-confidence available annotation track, but the attached code loops through every available listed track, appends every region set, and unions them. The README later also says tracks are merged, making the documentation internally inconsistent. The processed export contains **3,169 multi-track rows**, which is compatible with union behavior, although the actual production script is missing.

The original supplied script also cannot distinguish an API failure from “no annotation” in its returned interval list and does not emit `query_status` or `annotation_tracks`. The processed master handles missing results more carefully, but the code that did so was not supplied.

Recommendation: decide and document one method before using disorder fraction as a simulation QC variable:

1. highest-priority-track only, or
2. union of all listed tracks.

Then archive the exact script, package versions, query date, full result TSV, and ideally raw responses or a bulk-download version.

## 10. Identifier-mapping audit

Exact accession matching finds **1,264 of 1,770** McShane ED/NED records anywhere in the Jiang 20,588-protein universe; 506 do not match exactly. Of those 506, **399** have at least one accession from the Jiang universe inside McShane's grouped FASTA header. For 276 records, at least one such candidate also passes pLDDT ≥70.

Those candidates were **not** silently substituted. McShane FASTA headers can describe protein groups supported by shared peptides, so propagating one group-level decay label to a canonical accession or isoform requires a written mapping rule and review. Four records have more than one pLDDT-qualified candidate. All candidates are listed in `results/mcshane_grouped_header_candidate_mappings.csv`.

## 11. Manual source-to-output traces

Twenty-five deterministic traces were checked: four final proteins from each quadrant, one example for each major exclusion, one missing-MobiDB example per quadrant, and Q8NI36. Detailed CCBond counts and all derived fields are in `results/manual_source_to_output_traces.csv`.

| # | UniProt | pLDDT | Decay | Entanglement | Knot | FASTA | Membrane | MobiDB status | IDR fraction | Result |
| --- | --- | ---: | --- | --- | --- | --- | --- | --- | ---: | --- |
| 1 | A0AVT1 | 91.46 | ED | Native entangled | No | Yes | No | ok | 0.0466 | Included |
| 2 | A5YKK6 | 74.53 | ED | Native entangled | No | Yes | No | ok | 0.2597 | Included |
| 3 | A6NHR9 | 80.88 | ED | Native entangled | No | Yes | No | ok | 0.0464 | Included |
| 4 | O00170 | 91.09 | ED | Native entangled | No | Yes | No | ok | 0.1424 | Included |
| 5 | O00232 | 78.74 | ED | Not entangled | No | Yes | No | ok | 0.0417 | Included |
| 6 | O00410 | 92.03 | ED | Not entangled | No | Yes | No | ok | 0.0319 | Included |
| 7 | O00458 | 83.77 | ED | Not entangled | No | Yes | No | ok | 0.2550 | Included |
| 8 | O00505 | 86.45 | ED | Not entangled | No | Yes | No | ok | 0.1536 | Included |
| 9 | O00303 | 73.77 | NED | Native entangled | No | Yes | No | ok | 0.2465 | Included |
| 10 | O00487 | 81.85 | NED | Native entangled | No | Yes | No | ok | 0.0742 | Included |
| 11 | O43837 | 87.89 | NED | Native entangled | No | Yes | No | ok | 0.1221 | Included |
| 12 | O60264 | 74.94 | NED | Native entangled | No | Yes | No | ok | 0.5542 | Included |
| 13 | O15511 | 92.41 | NED | Not entangled | No | Yes | No | ok | 0.0927 | Included |
| 14 | O60783 | 86.82 | NED | Not entangled | No | Yes | No | ok | 0.2188 | Included |
| 15 | O75718 | 88.13 | NED | Not entangled | No | Yes | No | ok | 0.0998 | Included |
| 16 | P24844 | 83.25 | NED | Not entangled | No | Yes | No | ok | 0.1395 | Included |
| 17 | O75976 | 82.91 | ED | Covalent-lasso excluded | No | Yes | Yes | ok | 0.1587 | Excluded: covalent |
| 18 | O95905 | 76.28 | ED | Native entangled | Yes | Yes | No | ok | 0.2919 | Excluded: knot |
| 19 | P63092 | 91.28 | NED | Native entangled | No | No | Yes | ok | 0.1218 | Excluded: FASTA unavailable |
| 20 | O00116 | 88.35 | ED | Native entangled | No | Yes | Yes | ok | 0.1307 | Excluded: membrane |
| 21 | O00442 | 95.90 | ED | Native entangled | No | Yes | No | no annotation | — | Missing MobiDB |
| 22 | O00231 | 83.04 | ED | Not entangled | No | Yes | No | no annotation | — | Missing MobiDB |
| 23 | P06744 | 97.84 | NED | Native entangled | No | Yes | No | no annotation | — | Missing MobiDB |
| 24 | P39019 | 91.92 | NED | Not entangled | No | Yes | No | no annotation | — | Missing MobiDB |
| 25 | Q8NI36 | 83.73 | NED | Native entangled | No | Yes | No | ok | 0.0883 | Included; coordinate warning |

## 12. Exact discrepancies and ambiguity categories

| Category | Finding | Consequence |
| --- | --- | --- |
| Missing deliverables | Three generated Excel workbooks absent | Cannot verify exact workbook cells, sheets, ordering, or 606 workbook ID set. |
| Container integrity | Both repository ZIP copies invalid/truncated | Recovered source files are analyzable but not verifiable against an official archive hash. |
| Count/provenance | Reported 123 broad-master knots vs 122 from supplied source | Does not change the reproduced 606, but should be reconciled. |
| MobiDB provenance | Full TSV, production script, and raw responses absent | MobiDB values cannot be independently regenerated from supplied files. |
| Method definition | README says highest priority; code unions tracks | Disorder fractions depend on an unresolved methodological choice. |
| Coordinate versions | 18 interval/count display inconsistencies; Q8NI36 in 606 | Stored fractions follow clipping logic, but displayed coordinates and denominator are not version-compatible. |
| Identifier mapping | 506 ED/NED records lack exact Jiang accession; 399 have header candidates | Exact matching is reproducible but may lose valid coverage; automatic remapping is not justified. |
| Outcome retention | pLDDT-qualified UN records are blanked in broad master | Sensitivity analysis must re-merge McShane; it cannot use this column as-is. |
| Semantic wording | “Not entangled” means no detected NCLE record | Avoid interpreting it as definitive biological absence. |

## 13. Corrections and next actions

Before freezing the simulation table:

1. Attach/export the three exact Excel workbooks and rerun workbook comparison.
2. Archive the full production MobiDB TSV, exact production script, environment/package versions, and query date.
3. Decide whether MobiDB uses a highest-priority track or the union; make README and code agree, then regenerate if necessary.
4. Requery Q8NI36 using a version-matched sequence or retain it with a coordinate-mismatch flag until resolved.
5. Reconcile the 122-versus-123 knot count against the exact `AF4_knots.dat` version used on the original laptop.
6. Write and approve an accession/protein-group mapping policy before attempting to recover the 399 possible McShane-to-Jiang mappings.
7. Preserve raw `UN` labels in a separate master column for the planned sensitivity analysis.

## 14. Final conclusion on the 606 proteins

The **membership logic and counts of the 606-protein cohort are reproducible and internally coherent** from the available processed master and attached structural/QC sources. It can be used now as a provisional simulation-selection table, with Q8NI36 flagged.

For a final handoff or publication-quality archive, the 606 table is not fully verified until the exact workbook is attached, MobiDB provenance/method is frozen, and the one retained coordinate mismatch is resolved or explicitly documented.

## 15. Reproduction command

```bash
python audit/reproducibility_audit.py \
  --workspace /workspace/scratch/21898fbd418f \
  --output audit/results
```

The script is read-only with respect to supplied source/generated files and writes only new audit artifacts under `audit/results/`.
