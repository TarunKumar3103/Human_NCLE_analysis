# Protein decay input audit

This bundle accompanies `Protein_decay_project_review.md`. It audits the two supplied attachments without changing them or running molecular dynamics. The repository and experimental workbook are not duplicated in this bundle.

The tables describe the uploaded snapshot. They are not an approved matched simulation cohort.

| File | Contents |
|---|---|
| `audit_decay_inputs.py` | Recreates the protein-level aggregation, exact-leading-ID join, candidate mapping flags, and illustrative filter counts |
| `check_simulation_setup.py` | Checks a selected existing temperature-quench configuration and the local environment; submits no jobs |
| `audit/human_decay_join_audit.csv` | All 3,133 experimental entries with original fields and added curation flags |
| `audit/protein_level_structure_metadata.csv` | One row for each of 20,588 structural accessions |
| `audit/illustrative_filter_counts.csv` | Successive counts using selected old-repository filters |
| `audit/audit_summary.json` | Summary counts, audit policy, and SHA-256 hashes of the input datasets |
| `audit/local_simulation_preflight.json` | Results from the review environment, not Tarun's cluster environment |
| `audit/simulation_config_audit.json` | File-existence checks for the five setup references in each of 30 temperature-quench configurations |

To reproduce the input audit, use Python with pandas and openpyxl:

```bash
python audit_decay_inputs.py \
  --repo /path/to/Protein-entanglement-misfolding-determines-divergent-fates-main \
  --workbook '/path/to/1-s2.0-S009286741631248X-mmc4 (1).xlsx' \
  --out audit
```

To inspect one prepared example in the environment you intend to use for MD:

```bash
python check_simulation_setup.py \
  --repo /path/to/Protein-entanglement-misfolding-determines-divergent-fates-main \
  --group YU_E --protein P16152
```

These commands do not submit jobs or execute the original notebooks. The preflight locates optional simulation dependencies and attempts to list OpenMM platforms if OpenMM can be imported. The 30-configuration inventory was generated as a read-only review check; the per-protein preflight can independently repeat the setup-reference check for each listed protein.

Interpret the additional columns as follows:

- `source_excel_row`: original row in the `FinalTable_RPE_` worksheet; row 1 is the header.
- `accession`: original leading experimental accession. Alternative candidates never replace it automatically.
- `exact_leading_accession_match`: whether that exact identifier is present in the entanglement metadata.
- `entanglement_status`: `NCLE`, `no_detected_entanglement`, or `exclude_covalent_lasso` under the old notebook's aggregation convention. A missing joined value remains unknown.
- `n_entanglement_rows`, `n_noncovalent_rows`, `n_covalent_rows`: counts of supplied entanglement entries, not necessarily independent physical events.
- `sequence_length`: length in the supplied human FASTA. `structure_length_in_ent_file` is unavailable for no-entanglement placeholders.
- `sequence_length_matches_ent_file`: comparison when both lengths exist. A blank does not imply agreement; equal length alone does not establish sequence identity.
- `header_base_candidates_for_review`: accessions found in the visible FASTA header after base-accession normalization and restriction to available structural identifiers. This is not a validated mapping and may miss information truncated from the header.
- `binary_decay_for_review`: ED=0, NED=1, UN blank.
- `native_entangled_for_review`: no detected entanglement=0, NCLE=1; excluded or missing status blank.
- `passes_illustrative_basic_filters`: ED/NED, exact leading match, no covalent-lasso exclusion, no human knot-list entry, mean pLDDT at least 85, sequence length available, and absent from the supplied membrane list. It does not apply the additional length/asphericity steps and does not certify biological eligibility.
- `additional_review_required`: unresolved checks that still apply to the candidate pool, including secretory/disorder status and global matching measures.
- `distinct_ACO_values`: warns that the source ACO field can vary within a protein; no protein-wide ACO was manufactured.

In the first rows of the filter table, the four binary cell counts need not sum to the total: unresolved or excluded entanglement status and the UN decay class are retained at those stages.

Missing data were not imputed as zero, UN was not assigned to ED, and header candidates were not expanded into independent observations. The provisional cutoffs follow the older repository for illustration and need agreement for this study. No statistical significance test was run, and no production simulation configuration was created.
