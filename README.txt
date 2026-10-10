CURRENT LOCAL VERIFICATION
Read REPRODUCTION_RECORD.txt for completed checks and recovered sources.
SOURCE_MANIFEST_LOCAL.json records the current input-file hashes.

HUMAN NCLE-NED ANALYSIS: REPRODUCTION PACKAGE
Prepared for Tarun Kumar and Quyen Vu, 9 October 2026

START HERE
This package reproduces the reported human crude ORs, Model 1, and Model 2
from the saved protein-level analysis workbook. It also reconstructs the
saved complex-membership decisions, audits raw cohort filters through
membrane exclusion, and checks the processed MobiDB disorder table.
No MobiDB fraction or availability filter is used for statistical analyses.

QUICK START (Mac; commands go in VS Code's TERMINAL, not a notebook cell)
1. Extract the ZIP and open the Human_NCLE_reproduction folder in VS Code.
   The Terminal must show this folder; typing ls should show scripts and data.
2. Use Python 3.12 for the tested environment. Run each line in order:

   python3.12 -m venv .venv
   source .venv/bin/activate
   python -m pip install -r requirements.txt
   python scripts/export_complex_annotations.py
   python scripts/reproduce_results.py
   python scripts/audit_raw_sources.py
   python scripts/audit_disorder.py
   python scripts/prepare_mobidb_ids.py

If python3.12 is unavailable, use an installed compatible Python executable
to create the environment and record its version. The tested runtime is
Python 3.12.14. Dependency installation requires internet; the analysis and
audit commands above use local files and do not refresh online annotations.

EXPECTED CHECKS
- Recorded-decision reconstruction: PASS (1191 evidence records; 205 reviewed;
  28 provisional complex calls and 177 unresolved).
- Numerical reproduction: PASS. All 12 reported result rows match at numeric
  tolerance 1e-8; six Model 2 slope rows agree with an independent optimizer.
- Primary cohort 705; sensitivity cohort 1191; Model 2 samples 528 and 901.
- Source flow before secreted exclusion: primary 725, sensitivity 1230.
- Processed MobiDB master: 13690 rows, 12609 fractions available, 18 interval
  coordinate warnings. Primary statistics cohort still has 705 proteins;
  592 have a fraction. A missing fraction does not remove an analysis row.
- Average pLDDT >=70 query list: 13690 unique IDs.

DIRECTORY GUIDE
scripts/             Runnable current scripts: recovered originals and new
                     reproduction utilities, distinguished in the document.
data/                Unmodified reviewed reference workbook and processed
                     broad MobiDB export, plus exported saved complex calls.
raw_sources/         McShane human supplement, Jiang NCLE/pLDDT tables, knots,
                     proteome FASTA and historical membrane-keyword FASTA.
results/             Checked outputs from this package's validation run.
historical_audits/    Older source audits, repository notebook, Quyen's sample
                     disorder script, and simulation preflight. Historical
                     cutoffs/cohorts are not the current analysis protocol.
documentation/       Detailed manual, full text code collection, and an
                     independent-review prompt.
SOURCE_MANIFEST.json File hashes and provenance notes.

WHAT IS AND IS NOT RECOVERED
statistical_analysis.py is a recovered original, unchanged. The archived
production disorder downloader and Quyen's original sample are preserved.
The reproduction runner and audit utilities are newly written and tested.
The original code that converted UniProt text into initial complex classes
was not recovered. Saved UniProt evidence plus saved manual-review decisions
are supplied; export_complex_annotations.py reconstructs those recorded
binary labels exactly. It does not independently establish their biology.

The approved historical secreted-keyword source export was not recovered.
The final workbook preserves the selected cohorts and flags, so all reported
models reproduce. Full independent final-cohort construction from raw sources
needs that export. Do not derive a secreted list from eligible IDs: that would
make the validation circular. A current UniProt download is a new data version.

TO AUDIT THE FINAL SECRETED FILTER, ONCE THE HISTORICAL FILE IS AVAILABLE
The file must contain only secreted-positive accessions: FASTA, a JSON list,
or a one-ID-per-line text file. A CSV containing both 0 and 1 status rows
must first be restricted to positive rows. Pass its actual path:

   python scripts/audit_raw_sources.py --secreted-source raw_sources/secreted_ids.json --out results/full_source_audit

Expected final primary 705, sensitivity 1191, all flow and ID comparisons PASS.
There is intentionally no invented secreted_ids.json in this package.

The manual also supplies the complete original statistical_analysis.py command
for a raw-source rerun after recovering that source and exporting saved complex
annotations. Preserve the original reference workbook; write to a new folder.

OPTIONAL LIVE DISORDER QUERY
This is not required to reproduce association results. A live rerun can change
because UniProt/MobiDB sequences and annotations have changed. Its archived
script unions available tracks; the older README's priority wording should
not be interpreted as selecting one highest-priority track.

   python -m pip install -r requirements-live-query.txt
   python scripts/prepare_mobidb_ids.py
   python scripts/query_mobidb_disorder_production.py -f results/query_inputs/uniprots_test5.txt -o results/disorder_live_test5.tsv

Only after checking the five outputs, a full query uses:
   python scripts/query_mobidb_disorder_production.py -f results/query_inputs/uniprots_plddt_ge70.txt -o results/disorder_live_ge70.tsv

Do not overwrite archived inputs with refreshed outputs. Successful rows are
resumed; failed/no-annotation rows are retried. The recovered archived script
does not include a confirmed copy of any later laptop-specific IPv4 patch.

INTERPRETATION
The current cohorts do not show statistically detectable positive NCLE-NED
association in crude or adjusted analyses. This is not proof of no effect.
Model 2 is limited by missing complex annotations, provisional positive calls,
and one NED protein among explicit monomers. The complex OR is imprecise.
Mouse/yeast analyses and refolding MD are not completed by this package.

COPYABLE CODE
The Word manual contains the complete seven current scripts and Quyen's sample
script as selectable text. Individual .py files are the safest executable
copies. CODE_COMPENDIUM.txt includes those and all recovered historical Python
audit scripts; it is intended for reading/review, not execution as one file.
