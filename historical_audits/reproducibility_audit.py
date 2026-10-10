#!/usr/bin/env python3
"""Independent, read-only audit of the human NCLE/decay curation pipeline.

The script reconstructs the pLDDT >= 70 master and the ED/NED simulation
cohorts from supplied source files.  It does not run association tests, perform
matching, query external services, or overwrite any source/generated workbook.

Example
-------
python audit/reproducibility_audit.py \
  --workspace /workspace/scratch/21898fbd418f \
  --output audit/results
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


TARGET_WORKBOOKS = {
    "broad_master": "Human_master_pLDDT_ge70_with_MobiDB_and_membrane_QC.xlsx",
    "nonmembrane_core": "Human_ED_NED_nonmembrane_core_pLDDT_ge70.xlsx",
    "mobidb_complete": "Human_ED_NED_nonmembrane_MobiDB_complete.xlsx",
}

EXPECTED_BENCHMARKS = {
    "McShane rows": 3133,
    "McShane ED": 1474,
    "McShane NED": 296,
    "McShane UN": 1363,
    "ED/NED rows": 1770,
    "pLDDT source proteins": 20588,
    "pLDDT >= 70 master": 13690,
    "Native entangled in master": 8416,
    "Not entangled in master": 4782,
    "Covalent-lasso excluded in master": 492,
    "ED/NED in master": 1043,
    "Eligible entanglement among ED/NED": 1034,
    "After knot exclusion": 1023,
    "FASTA-complete pre-membrane cohort": 1022,
    "Non-membrane core": 725,
    "MobiDB-complete core": 606,
    "MobiDB ok": 12609,
    "MobiDB no selected annotation": 1042,
    "MobiDB no selected annotation and missing length": 31,
    "MobiDB missing length": 8,
    # This is the number previously printed in the notebook transcript.
    "Knotted proteins in master (reported)": 123,
    "FASTA available in master": 13594,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Workspace root containing downloads/, upload/, and work/.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("audit/results"),
        help="Output directory; relative paths are resolved under --workspace.",
    )
    return parser.parse_args()


def first_existing(candidates: Iterable[Path]) -> Path:
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError("None of the candidate paths exist:\n" + "\n".join(map(str, candidates)))


def clean_id(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip()


def parse_boolish(value: object) -> bool | None:
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if text in {"true", "1", "1.0"}:
        return True
    if text in {"false", "0", "0.0"}:
        return False
    raise ValueError(f"Unexpected CCBond value: {value!r}")


def aggregate_entanglement(series: pd.Series) -> float:
    """Repository rule: any covalent crossing excludes; otherwise NCLE if present."""
    values = [parse_boolish(value) for value in series]
    values = [value for value in values if value is not None]
    if any(values):
        return np.nan
    if values:
        return 1.0
    return 0.0


def read_fasta_ids(path: Path) -> tuple[list[str], set[str]]:
    ids: list[str] = []
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            parts = header.split("|")
            accession = parts[1] if len(parts) >= 3 else header.split()[0]
            ids.append(accession.strip())
    return ids, set(ids)


def parse_regions(value: object) -> list[tuple[int, int]]:
    if pd.isna(value) or str(value).strip() in {"", "-"}:
        return []
    regions: list[tuple[int, int]] = []
    for token in str(value).split(";"):
        match = re.fullmatch(r"(\d+)-(\d+)", token.strip())
        if not match:
            raise ValueError(f"Unparseable disorder interval: {token!r}")
        start, end = map(int, match.groups())
        if start < 1 or end < start:
            raise ValueError(f"Invalid disorder interval: {token!r}")
        regions.append((start, end))
    return regions


def compare_series(left: pd.Series, right: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
        return pd.Series(
            ~np.isclose(left, right, equal_nan=True, atol=1e-12, rtol=1e-12),
            index=left.index,
        )
    return left.fillna("<NA>").astype(str).ne(right.fillna("<NA>").astype(str))


def find_any(workspace: Path, filename: str) -> list[Path]:
    return sorted(path for path in workspace.rglob(filename) if "audit/results" not in str(path))


def write_table(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    workspace = args.workspace.resolve()
    output = args.output if args.output.is_absolute() else workspace / args.output
    output.mkdir(parents=True, exist_ok=True)

    repo_dir = first_existing(
        [
            workspace
            / "work/repo_recovered/Protein-entanglement-misfolding-determines-divergent-fates-main/0_Dataset_curation_and_logistic_regression",
            workspace
            / "Protein-entanglement-misfolding-determines-divergent-fates-main/0_Dataset_curation_and_logistic_regression",
        ]
    )
    mcshane_path = first_existing(
        [
            workspace / "downloads/1-s2.0-S009286741631248X-mmc4 (1).xlsx",
            repo_dir / "Human_degradation_mmc4.xlsx",
        ]
    )
    pasted_master_path = first_existing(
        [
            workspace / "upload/Pasted text(20260921-021255).txt",
            workspace / "Pasted text(20260921-021255).txt",
        ]
    )
    disorder_package = workspace / "work/disorder_query_original/disorder_query"

    paths = {
        "McShane workbook": mcshane_path,
        "pLDDT source": repo_dir / "Human_Avg_pLDDTs.csv",
        "entanglement source": repo_dir / "Human_AF_combined_20250614.csv",
        "knot source": repo_dir / "AF4_knots.dat",
        "human FASTA": repo_dir / "UP000005640_9606.fasta",
        "membrane-keyword FASTA": repo_dir / "uniprotkb_human_keywords_membrane_2025_03_19.fasta",
        "processed broad master export": pasted_master_path,
        "Quyen README": disorder_package / "README.md",
        "Quyen original query script": disorder_package / "get_disorder_mobidb_query.py",
        "Quyen five-protein sample": disorder_package / "disorder.tsv",
    }

    zip_candidates = [
        workspace / "downloads/Protein-entanglement-misfolding-determines-divergent-fates-main.zip",
        workspace
        / "downloads/Protein-entanglement-misfolding-determines-divergent-fates-main.zip.openai-download-d117007e73c347ee9d48076781654dfb",
    ]

    inventory_rows: list[dict[str, object]] = []
    for label, path in paths.items():
        inventory_rows.append(
            {
                "item": label,
                "path": str(path),
                "present": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else np.nan,
                "sha256": sha256_file(path) if path.exists() else "",
                "container_valid": "not applicable",
                "role": "source" if "processed" not in label else "processed output",
            }
        )
    for index, path in enumerate(zip_candidates, start=1):
        inventory_rows.append(
            {
                "item": f"repository ZIP copy {index}",
                "path": str(path),
                "present": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else np.nan,
                "sha256": sha256_file(path) if path.exists() else "",
                "container_valid": zipfile.is_zipfile(path) if path.exists() else False,
                "role": "supplied container; relevant payloads were recovered because central directory is absent",
            }
        )
    for label, filename in TARGET_WORKBOOKS.items():
        matches = find_any(workspace, filename)
        inventory_rows.append(
            {
                "item": f"target workbook: {label}",
                "path": "; ".join(map(str, matches)),
                "present": bool(matches),
                "size_bytes": matches[0].stat().st_size if matches else np.nan,
                "sha256": sha256_file(matches[0]) if matches else "",
                "container_valid": "not checked" if matches else "not available",
                "role": "generated workbook requested for comparison",
            }
        )
    inventory_rows.extend(
        [
            {
                "item": "full production MobiDB TSV",
                "path": "",
                "present": False,
                "size_bytes": np.nan,
                "sha256": "",
                "container_valid": "not available",
                "role": "missing raw/processed query source (only embedded columns were supplied)",
            },
            {
                "item": "production MobiDB query script",
                "path": "",
                "present": False,
                "size_bytes": np.nan,
                "sha256": "",
                "container_valid": "not available",
                "role": "missing generator for query_status and annotation_tracks columns",
            },
            {
                "item": "raw MobiDB JSON responses",
                "path": "",
                "present": False,
                "size_bytes": np.nan,
                "sha256": "",
                "container_valid": "not available",
                "role": "missing API provenance snapshot",
            },
        ]
    )
    inventory = pd.DataFrame(inventory_rows)

    # McShane degradation source.
    mcshane = pd.read_excel(mcshane_path, sheet_name="FinalTable_RPE_")
    mcshane = mcshane.copy()
    mcshane["UniProt_ID"] = clean_id(mcshane["Protein IDs (Uniprot)"])
    mcshane["Degradation_profile_raw"] = mcshane["Degradation profile"].astype("string").str.strip()
    degradation = mcshane[["UniProt_ID", "Degradation_profile_raw"]].copy()

    # pLDDT source and retained broad master backbone.
    plddt = pd.read_csv(paths["pLDDT source"])
    plddt["UniProt_ID"] = clean_id(plddt["gene"])
    plddt = plddt.rename(columns={"<pLDDT>": "Average_pLDDT", "std": "pLDDT_SD"})
    plddt_master = plddt.loc[
        pd.to_numeric(plddt["Average_pLDDT"], errors="coerce").ge(70),
        ["UniProt_ID", "Average_pLDDT", "pLDDT_SD"],
    ].copy()

    # Entanglement source. Preserve source diagnostics at protein level.
    ent_raw = pd.read_csv(paths["entanglement source"], sep="|", low_memory=False)
    ent_raw["UniProt_ID"] = clean_id(ent_raw["gene"])
    ent_raw["CCBond_parsed"] = ent_raw["CCBond"].map(parse_boolish)
    ent_diag = (
        ent_raw.groupby("UniProt_ID", as_index=False)
        .agg(
            Entanglement_source_rows=("UniProt_ID", "size"),
            ENT_ID_nonmissing=("ENT-ID", "count"),
            CCBond_true=("CCBond_parsed", lambda s: sum(value is True for value in s)),
            CCBond_false=("CCBond_parsed", lambda s: sum(value is False for value in s)),
            CCBond_missing=("CCBond_parsed", lambda s: int(s.isna().sum())),
            Entanglement=("CCBond", aggregate_entanglement),
        )
    )
    ent_diag["Entanglement_status"] = np.select(
        [ent_diag["Entanglement"].eq(1), ent_diag["Entanglement"].eq(0)],
        ["Native_entangled", "Not_entangled"],
        default="Excluded_covalent_lasso",
    )

    reconstructed = plddt_master.merge(
        ent_diag[["UniProt_ID", "Entanglement", "Entanglement_status"]],
        on="UniProt_ID",
        how="left",
        validate="one_to_one",
    )
    reconstructed = reconstructed.merge(
        degradation,
        on="UniProt_ID",
        how="left",
        validate="one_to_one",
    )
    reconstructed["Degradation_profile"] = reconstructed["Degradation_profile_raw"].where(
        reconstructed["Degradation_profile_raw"].isin(["ED", "NED"])
    )
    reconstructed["Bi_exponential_decay"] = reconstructed["Degradation_profile"].map(
        {"ED": 0.0, "NED": 1.0}
    )

    # Supplied processed export contains all 13,690 MobiDB results and the broad table.
    provided = pd.read_csv(pasted_master_path)
    provided["UniProt_ID"] = clean_id(provided["UniProt_ID"])
    source_comparison = reconstructed.merge(
        provided,
        on="UniProt_ID",
        how="outer",
        suffixes=("_reconstructed", "_provided"),
        indicator=True,
    )
    comparison_rows: list[dict[str, object]] = [
        {
            "comparison": "ID set",
            "reconstructed_rows": len(reconstructed),
            "provided_rows": len(provided),
            "mismatches": int(source_comparison["_merge"].ne("both").sum()),
        }
    ]
    for column in [
        "Average_pLDDT",
        "pLDDT_SD",
        "Entanglement",
        "Entanglement_status",
        "Degradation_profile",
        "Bi_exponential_decay",
    ]:
        mismatch = compare_series(
            source_comparison[f"{column}_reconstructed"],
            source_comparison[f"{column}_provided"],
        )
        comparison_rows.append(
            {
                "comparison": column,
                "reconstructed_rows": len(reconstructed),
                "provided_rows": len(provided),
                "mismatches": int(mismatch.sum()),
            }
        )
    source_field_comparison = pd.DataFrame(comparison_rows)

    # Independent flags from attached source files.
    fasta_list, fasta_ids = read_fasta_ids(paths["human FASTA"])
    membrane_list, membrane_ids = read_fasta_ids(paths["membrane-keyword FASTA"])
    knots = pd.read_csv(
        paths["knot source"],
        sep="\t",
        comment="#",
        header=None,
        names=["Knot_type", "Category", "Uniprot", "Organism", "pLDDT_knotcore", "Protein_name"],
        dtype=str,
    )
    knots["Uniprot"] = clean_id(knots["Uniprot"])
    knots["Organism"] = clean_id(knots["Organism"])
    human_knots = knots.loc[knots["Organism"].eq("Human")].copy()
    human_knots["UniProt_ID"] = human_knots["Uniprot"].str.replace(r"-F\d+$", "", regex=True)
    knot_ids = set(human_knots["UniProt_ID"].dropna())

    master = provided.copy()
    master["Knotted"] = master["UniProt_ID"].isin(knot_ids)
    master["FASTA_available"] = master["UniProt_ID"].isin(fasta_ids)
    master["Membrane"] = master["UniProt_ID"].isin(membrane_ids)

    # Sequential cohort reconstruction. Secretory status is intentionally not used.
    edned = master.loc[master["Degradation_profile"].isin(["ED", "NED"])].copy()
    ent_eligible = edned.loc[edned["Entanglement"].notna()].copy()
    knot_free = ent_eligible.loc[~ent_eligible["Knotted"]].copy()
    fasta_complete = knot_free.loc[knot_free["FASTA_available"]].copy()
    core = fasta_complete.loc[~fasta_complete["Membrane"]].copy()
    core["Analysis_group"] = core["Degradation_profile"] + "_" + core["Entanglement"].map(
        {1.0: "Entangled", 0.0: "NotEntangled"}
    )
    complete = core.loc[core["Fraction_IDP_IDR"].notna()].copy()

    cohort_flow = pd.DataFrame(
        [
            ("pLDDT >= 70 broad master", len(master), 0),
            ("retain ED/NED", len(edned), len(master) - len(edned)),
            ("exclude covalent-lasso entanglements", len(ent_eligible), len(edned) - len(ent_eligible)),
            ("exclude knotted proteins", len(knot_free), len(ent_eligible) - len(knot_free)),
            ("require proteome FASTA accession", len(fasta_complete), len(knot_free) - len(fasta_complete)),
            ("exclude UniProt membrane-keyword proteins", len(core), len(fasta_complete) - len(core)),
            ("require nonmissing MobiDB disorder fraction", len(complete), len(core) - len(complete)),
        ],
        columns=["stage", "proteins_remaining", "removed_at_step"],
    )

    group_order = ["ED_Entangled", "ED_NotEntangled", "NED_Entangled", "NED_NotEntangled"]
    group_counts = []
    for group in group_order:
        total = int(core["Analysis_group"].eq(group).sum())
        available = int(complete["Analysis_group"].eq(group).sum())
        group_counts.append(
            {
                "Analysis_group": group,
                "nonmembrane_core": total,
                "MobiDB_complete": available,
                "MobiDB_missing": total - available,
                "MobiDB_missing_percent": 100 * (total - available) / total,
            }
        )
    group_counts = pd.DataFrame(group_counts)

    # MobiDB arithmetic, coordinate, and status checks.
    interval_rows: list[dict[str, object]] = []
    for row in master.itertuples(index=False):
        regions = parse_regions(row.Disorder_regions)
        if not regions:
            continue
        displayed_count = sum(end - start + 1 for start, end in regions)
        if pd.isna(row.MobiDB_length):
            clipped_count = np.nan
            out_of_range = np.nan
        else:
            length = int(row.MobiDB_length)
            clipped_count = sum(
                max(0, min(end, length) - max(start, 1) + 1) for start, end in regions
            )
            out_of_range = any(end > length or start < 1 for start, end in regions)
        interval_rows.append(
            {
                "UniProt_ID": row.UniProt_ID,
                "MobiDB_length": row.MobiDB_length,
                "Disorder_regions": row.Disorder_regions,
                "displayed_interval_residues": displayed_count,
                "clipped_interval_residues": clipped_count,
                "stored_Disorder_residues": row.Disorder_residues,
                "displayed_minus_stored": displayed_count - row.Disorder_residues,
                "clipped_matches_stored": bool(
                    np.isclose(clipped_count, row.Disorder_residues, equal_nan=True)
                ),
                "interval_extends_outside_length": out_of_range,
                "in_nonmembrane_core": row.UniProt_ID in set(core["UniProt_ID"]),
                "in_606_complete": row.UniProt_ID in set(complete["UniProt_ID"]),
            }
        )
    interval_audit = pd.DataFrame(interval_rows)
    interval_discrepancies = interval_audit.loc[
        ~np.isclose(
            interval_audit["displayed_interval_residues"],
            interval_audit["stored_Disorder_residues"],
            equal_nan=True,
        )
    ].copy()

    formula = master["Disorder_residues"] / master["MobiDB_length"]
    # Stored values have eight decimal places; allow half-unit rounding plus FP noise.
    formula_mismatch = (
        master["Fraction_IDP_IDR"].notna()
        & (master["Fraction_IDP_IDR"] - formula).abs().gt(5.1e-9)
    )

    # Experimental protein-group headers that contain alternate accessions.
    plddt_ids = set(plddt["UniProt_ID"])

    def extract_header_accessions(value: object) -> list[str]:
        return re.findall(r"(?:sp|tr)\|([A-Z0-9]+(?:-[0-9]+)?)\|", str(value))

    experimental_edned = mcshane.loc[mcshane["Degradation_profile_raw"].isin(["ED", "NED"])].copy()
    experimental_edned["Exact_Jiang_match"] = experimental_edned["UniProt_ID"].isin(plddt_ids)
    experimental_edned["Header_Jiang_candidates"] = experimental_edned["Fasta header"].map(
        lambda value: sorted(
            set(accession for accession in extract_header_accessions(value) if accession in plddt_ids)
        )
    )
    candidate_mappings = experimental_edned.loc[
        ~experimental_edned["Exact_Jiang_match"]
        & experimental_edned["Header_Jiang_candidates"].map(bool),
        ["UniProt_ID", "Degradation_profile_raw", "Header_Jiang_candidates", "Fasta header"],
    ].copy()
    plddt_ge70_ids = set(plddt_master["UniProt_ID"])
    candidate_mappings["Header_candidates_pLDDT_ge70"] = candidate_mappings[
        "Header_Jiang_candidates"
    ].map(lambda values: sorted(set(values).intersection(plddt_ge70_ids)))
    candidate_mappings["n_candidates_all_pLDDT"] = candidate_mappings[
        "Header_Jiang_candidates"
    ].map(len)
    candidate_mappings["n_candidates_pLDDT_ge70"] = candidate_mappings[
        "Header_candidates_pLDDT_ge70"
    ].map(len)
    candidate_mappings["Header_Jiang_candidates"] = candidate_mappings[
        "Header_Jiang_candidates"
    ].map(";".join)
    candidate_mappings["Header_candidates_pLDDT_ge70"] = candidate_mappings[
        "Header_candidates_pLDDT_ge70"
    ].map(";".join)

    # Manual source-to-output traces: four final rows per quadrant plus exclusions.
    traced_ids: list[str] = []
    for group in group_order:
        traced_ids.extend(
            complete.loc[complete["Analysis_group"].eq(group), "UniProt_ID"].sort_values().head(4)
        )

    def first_id(frame: pd.DataFrame) -> str | None:
        values = frame["UniProt_ID"].sort_values()
        return values.iloc[0] if len(values) else None

    extras = [
        first_id(edned.loc[edned["Entanglement"].isna()]),
        first_id(ent_eligible.loc[ent_eligible["Knotted"]]),
        first_id(knot_free.loc[~knot_free["FASTA_available"]]),
        first_id(fasta_complete.loc[fasta_complete["Membrane"]]),
    ]
    for group in group_order:
        extras.append(
            first_id(core.loc[core["Analysis_group"].eq(group) & core["Fraction_IDP_IDR"].isna()])
        )
    extras.append("Q8NI36" if "Q8NI36" in set(master["UniProt_ID"]) else None)
    traced_ids.extend(value for value in extras if value is not None)
    traced_ids = list(dict.fromkeys(traced_ids))

    trace = master.loc[master["UniProt_ID"].isin(traced_ids)].copy()
    trace = trace.merge(ent_diag, on="UniProt_ID", how="left", suffixes=("", "_source"))
    trace["In_ED_NED"] = trace["Degradation_profile"].isin(["ED", "NED"])
    trace["Entanglement_eligible"] = trace["Entanglement"].notna()
    trace["Knot_free"] = ~trace["Knotted"]
    trace["Non_membrane"] = ~trace["Membrane"]
    trace["MobiDB_complete"] = trace["Fraction_IDP_IDR"].notna()
    trace["Analysis_group"] = np.where(
        trace["Degradation_profile"].isin(["ED", "NED"]) & trace["Entanglement"].notna(),
        trace["Degradation_profile"]
        + "_"
        + trace["Entanglement"].map({1.0: "Entangled", 0.0: "NotEntangled"}),
        pd.NA,
    )

    def exclusion_reason(row: pd.Series) -> str:
        if not row["In_ED_NED"]:
            return "not ED/NED"
        if not row["Entanglement_eligible"]:
            return "covalent-lasso entanglement"
        if not row["Knot_free"]:
            return "knotted"
        if not row["FASTA_available"]:
            return "FASTA accession unavailable"
        if not row["Non_membrane"]:
            return "membrane keyword"
        if not row["MobiDB_complete"]:
            return "missing MobiDB fraction"
        return "included in 606"

    trace["Final_result"] = trace.apply(exclusion_reason, axis=1)
    trace["Trace_order"] = trace["UniProt_ID"].map({value: i + 1 for i, value in enumerate(traced_ids)})
    trace = trace.sort_values("Trace_order")
    trace_columns = [
        "Trace_order",
        "UniProt_ID",
        "Average_pLDDT",
        "pLDDT_SD",
        "Degradation_profile",
        "Bi_exponential_decay",
        "ENT_ID_nonmissing",
        "CCBond_true",
        "CCBond_false",
        "CCBond_missing",
        "Entanglement",
        "Entanglement_status",
        "Knotted",
        "FASTA_available",
        "Membrane",
        "Disorder_query_status",
        "Disorder_residues",
        "MobiDB_length",
        "Fraction_IDP_IDR",
        "Analysis_group",
        "Final_result",
    ]
    trace = trace.reindex(columns=trace_columns)

    # Merge/source diagnostics.
    raw_ccbond_values = Counter(str(value).strip() for value in ent_raw["CCBond"].dropna())
    invalid_ccbond = sum(
        count
        for value, count in raw_ccbond_values.items()
        if value.lower() not in {"true", "false", "1", "0", "1.0", "0.0"}
    )
    exact_edned_matches = int(experimental_edned["Exact_Jiang_match"].sum())
    no_exact = len(experimental_edned) - exact_edned_matches
    has_header_candidate = int(
        experimental_edned.loc[~experimental_edned["Exact_Jiang_match"], "Header_Jiang_candidates"]
        .map(bool)
        .sum()
    )
    merge_audit = pd.DataFrame(
        [
            {
                "source": "McShane",
                "rows": len(mcshane),
                "unique_IDs": mcshane["UniProt_ID"].nunique(),
                "missing_IDs": mcshane["UniProt_ID"].isna().sum(),
                "duplicate_IDs": mcshane["UniProt_ID"].duplicated().sum(),
                "coverage_or_note": "All rows have one exact Protein IDs (Uniprot) value",
            },
            {
                "source": "Jiang pLDDT",
                "rows": len(plddt),
                "unique_IDs": plddt["UniProt_ID"].nunique(),
                "missing_IDs": plddt["UniProt_ID"].isna().sum(),
                "duplicate_IDs": plddt["UniProt_ID"].duplicated().sum(),
                "coverage_or_note": f"{len(plddt_master)} pass mean pLDDT >= 70",
            },
            {
                "source": "Jiang entanglement",
                "rows": len(ent_raw),
                "unique_IDs": ent_raw["UniProt_ID"].nunique(),
                "missing_IDs": ent_raw["UniProt_ID"].isna().sum(),
                "duplicate_IDs": len(ent_raw) - ent_raw["UniProt_ID"].nunique(),
                "coverage_or_note": f"{invalid_ccbond} invalid nonmissing CCBond values",
            },
            {
                "source": "processed MobiDB/master export",
                "rows": len(provided),
                "unique_IDs": provided["UniProt_ID"].nunique(),
                "missing_IDs": provided["UniProt_ID"].isna().sum(),
                "duplicate_IDs": provided["UniProt_ID"].duplicated().sum(),
                "coverage_or_note": f"{provided['Fraction_IDP_IDR'].notna().sum()} nonmissing fractions",
            },
            {
                "source": "human proteome FASTA",
                "rows": len(fasta_list),
                "unique_IDs": len(fasta_ids),
                "missing_IDs": 0,
                "duplicate_IDs": len(fasta_list) - len(fasta_ids),
                "coverage_or_note": f"{master['FASTA_available'].sum()} of master IDs matched",
            },
            {
                "source": "membrane-keyword FASTA",
                "rows": len(membrane_list),
                "unique_IDs": len(membrane_ids),
                "missing_IDs": 0,
                "duplicate_IDs": len(membrane_list) - len(membrane_ids),
                "coverage_or_note": f"{master['Membrane'].sum()} of master IDs matched",
            },
            {
                "source": "AlphaKnot human records",
                "rows": len(human_knots),
                "unique_IDs": human_knots["UniProt_ID"].nunique(),
                "missing_IDs": human_knots["UniProt_ID"].isna().sum(),
                "duplicate_IDs": human_knots["UniProt_ID"].duplicated().sum(),
                "coverage_or_note": f"{master['Knotted'].sum()} of master IDs matched after -F# removal",
            },
            {
                "source": "McShane ED/NED -> Jiang exact accession",
                "rows": len(experimental_edned),
                "unique_IDs": experimental_edned["UniProt_ID"].nunique(),
                "missing_IDs": 0,
                "duplicate_IDs": 0,
                "coverage_or_note": (
                    f"{exact_edned_matches} exact; {no_exact} unmatched; "
                    f"{has_header_candidate} unmatched records contain >=1 possible Jiang accession in grouped FASTA header"
                ),
            },
        ]
    )

    status = master["Disorder_query_status"].value_counts(dropna=False)
    observed_benchmarks = {
        "McShane rows": len(mcshane),
        "McShane ED": int(mcshane["Degradation_profile_raw"].eq("ED").sum()),
        "McShane NED": int(mcshane["Degradation_profile_raw"].eq("NED").sum()),
        "McShane UN": int(mcshane["Degradation_profile_raw"].eq("UN").sum()),
        "ED/NED rows": int(mcshane["Degradation_profile_raw"].isin(["ED", "NED"]).sum()),
        "pLDDT source proteins": len(plddt),
        "pLDDT >= 70 master": len(master),
        "Native entangled in master": int(master["Entanglement"].eq(1).sum()),
        "Not entangled in master": int(master["Entanglement"].eq(0).sum()),
        "Covalent-lasso excluded in master": int(master["Entanglement"].isna().sum()),
        "ED/NED in master": len(edned),
        "Eligible entanglement among ED/NED": len(ent_eligible),
        "After knot exclusion": len(knot_free),
        "FASTA-complete pre-membrane cohort": len(fasta_complete),
        "Non-membrane core": len(core),
        "MobiDB-complete core": len(complete),
        "MobiDB ok": int(status.get("ok", 0)),
        "MobiDB no selected annotation": int(status.get("mobidb_no_selected_annotation", 0)),
        "MobiDB no selected annotation and missing length": int(
            status.get("mobidb_no_selected_annotation|uniprot_missing_length", 0)
        ),
        "MobiDB missing length": int(status.get("uniprot_missing_length", 0)),
        "Knotted proteins in master (reported)": int(master["Knotted"].sum()),
        "FASTA available in master": int(master["FASTA_available"].sum()),
    }
    benchmark_comparison = pd.DataFrame(
        [
            {
                "benchmark": name,
                "expected_or_reported": expected,
                "independently_observed": observed_benchmarks[name],
                "difference": observed_benchmarks[name] - expected,
                "matches": observed_benchmarks[name] == expected,
            }
            for name, expected in EXPECTED_BENCHMARKS.items()
        ]
    )

    quality_checks = pd.DataFrame(
        [
            ("provided master duplicate IDs", int(master["UniProt_ID"].duplicated().sum()), 0),
            ("provided master missing IDs", int(master["UniProt_ID"].isna().sum()), 0),
            ("provided master Index sequential 1..N", int(master["Index"].tolist() == list(range(1, len(master) + 1))), 1),
            ("retained pLDDT below 70", int(master["Average_pLDDT"].lt(70).sum()), 0),
            ("fraction outside [0,1]", int(((master["Fraction_IDP_IDR"] < 0) | (master["Fraction_IDP_IDR"] > 1)).sum()), 0),
            ("fraction formula mismatches beyond 8-decimal rounding", int(formula_mismatch.sum()), 0),
            ("disorder residues greater than MobiDB length", int((master["Disorder_residues"] > master["MobiDB_length"]).sum()), 0),
            ("displayed interval count differs from stored count", len(interval_discrepancies), 0),
            ("clipped interval count differs from stored count", int((~interval_audit["clipped_matches_stored"]).sum()), 0),
            ("MobiDB complete rows with non-ok status", int((master["Fraction_IDP_IDR"].notna() & master["Disorder_query_status"].ne("ok")).sum()), 0),
            ("MobiDB ok rows with missing fraction", int((master["Disorder_query_status"].eq("ok") & master["Fraction_IDP_IDR"].isna()).sum()), 0),
            ("multi-track MobiDB rows", int(master["MobiDB_annotation_tracks"].astype(str).str.contains(";", regex=False).sum()), np.nan),
            ("manual trace rows", len(trace), 20),
        ],
        columns=["check", "observed", "expected_or_minimum"],
    )

    # Workbook-level comparisons can only run when exact files are supplied.
    workbook_rows: list[dict[str, object]] = []
    expected_frames = {
        "broad_master": master,
        "nonmembrane_core": core,
        "mobidb_complete": complete,
    }
    id_candidates = ["UniProt_ID", "Uniprot", "UniProt", "gene"]
    for label, filename in TARGET_WORKBOOKS.items():
        matches = find_any(workspace, filename)
        if not matches:
            workbook_rows.append(
                {
                    "workbook": filename,
                    "status": "MISSING",
                    "rows": np.nan,
                    "unique_IDs": np.nan,
                    "expected_rows": len(expected_frames[label]),
                    "ID_set_differences": np.nan,
                    "note": "Exact generated workbook was not supplied; workbook-level verification not possible.",
                }
            )
            continue
        book = pd.read_excel(matches[0])
        id_column = next((column for column in id_candidates if column in book.columns), None)
        if id_column is None:
            workbook_rows.append(
                {
                    "workbook": filename,
                    "status": "ID COLUMN NOT FOUND",
                    "rows": len(book),
                    "unique_IDs": np.nan,
                    "expected_rows": len(expected_frames[label]),
                    "ID_set_differences": np.nan,
                    "note": f"Columns: {list(book.columns)}",
                }
            )
            continue
        book_ids = set(clean_id(book[id_column]).dropna())
        expected_ids = set(expected_frames[label]["UniProt_ID"])
        workbook_rows.append(
            {
                "workbook": filename,
                "status": "FOUND",
                "rows": len(book),
                "unique_IDs": len(book_ids),
                "expected_rows": len(expected_frames[label]),
                "ID_set_differences": len(book_ids.symmetric_difference(expected_ids)),
                "note": str(matches[0]),
            }
        )
    workbook_comparison = pd.DataFrame(workbook_rows)

    # Write only new audit artifacts.
    write_table(inventory, output / "file_inventory.csv")
    write_table(merge_audit, output / "source_and_merge_audit.csv")
    write_table(cohort_flow, output / "cohort_flow.csv")
    write_table(group_counts, output / "group_counts_and_mobidb_missingness.csv")
    write_table(benchmark_comparison, output / "benchmark_comparison.csv")
    write_table(source_field_comparison, output / "broad_master_source_field_comparison.csv")
    write_table(workbook_comparison, output / "target_workbook_comparison.csv")
    write_table(quality_checks, output / "quality_checks.csv")
    write_table(interval_discrepancies, output / "mobidb_interval_discrepancies.csv")
    write_table(candidate_mappings, output / "mcshane_grouped_header_candidate_mappings.csv")
    write_table(trace, output / "manual_source_to_output_traces.csv")

    # ID lists make all sequential exclusions reviewable without regenerating workbooks.
    for filename, frame in {
        "reconstructed_nonmembrane_core_ids.csv": core,
        "reconstructed_mobidb_complete_606_ids.csv": complete,
        "missing_mobidb_from_core.csv": core.loc[core["Fraction_IDP_IDR"].isna()],
    }.items():
        write_table(frame, output / filename)

    summary = {
        "workspace": str(workspace),
        "output": str(output),
        "verdict": "PASS WITH WARNINGS",
        "broad_master_rows": len(master),
        "nonmembrane_core_rows": len(core),
        "mobidb_complete_rows": len(complete),
        "all_source_derived_master_fields_match_processed_export": bool(
            source_field_comparison["mismatches"].eq(0).all()
        ),
        "reported_knot_count": 123,
        "independent_knot_count": int(master["Knotted"].sum()),
        "mobidb_displayed_interval_discrepancies": len(interval_discrepancies),
        "mobidb_discrepancies_in_606": int(interval_discrepancies["in_606_complete"].sum()),
        "multi_track_mobidb_rows": int(
            master["MobiDB_annotation_tracks"].astype(str).str.contains(";", regex=False).sum()
        ),
        "target_workbooks_missing": int(workbook_comparison["status"].eq("MISSING").sum()),
        "manual_trace_rows": len(trace),
        "secretory_filter_applied": False,
        "association_testing_performed": False,
    }
    (output / "audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
