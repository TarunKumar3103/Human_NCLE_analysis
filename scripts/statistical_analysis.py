#!/usr/bin/env python3
"""Reconstruct and analyze the human NCLE/degradation cohorts.

Primary analysis: NED versus ED.
Sensitivity analysis: NED versus ED + UN.

MobiDB disorder is intentionally not used as an inclusion criterion.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mcshane", type=Path, required=True)
    parser.add_argument("--entanglement", type=Path, required=True)
    parser.add_argument("--plddt", type=Path, required=True)
    parser.add_argument("--knots", type=Path, required=True)
    parser.add_argument("--proteome-fasta", type=Path, required=True)
    parser.add_argument("--membrane-fasta", type=Path, required=True)
    parser.add_argument("--secreted-source", type=Path, required=True)
    parser.add_argument("--monomer-fasta", type=Path, action="append", default=[])
    parser.add_argument("--complex-annotation", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def normalize_id(value: object) -> str:
    return str(value).strip()


def merge_ccbond(series: pd.Series) -> float:
    values: list[bool] = []
    for item in series.dropna():
        if isinstance(item, (bool, np.bool_)):
            values.append(bool(item))
        elif str(item).strip() in {"True", "1", "1.0"}:
            values.append(True)
        elif str(item).strip() in {"False", "0", "0.0"}:
            values.append(False)
    if any(values):
        return np.nan
    if values:
        return 1.0
    return 0.0


def parse_fasta(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current: str | None = None
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                parts = line[1:].split("|")
                current = parts[1] if len(parts) >= 3 else parts[0].split()[0]
                current = normalize_id(current)
                lengths.setdefault(current, 0)
            elif current is not None:
                lengths[current] += len(line)
    return lengths


def parse_id_source(path: Path) -> set[str]:
    suffix = path.suffix.lower()
    if suffix in {".fa", ".faa", ".fasta"}:
        return set(parse_fasta(path))
    if suffix == ".json":
        payload = json.loads(path.read_text())
        values = payload.get("secreted", payload) if isinstance(payload, dict) else payload
        return {normalize_id(value) for value in values}
    values: set[str] = set()
    with path.open() as handle:
        for raw_line in handle:
            value = raw_line.strip().split("\t")[0].split(",")[0]
            if value and value.lower() not in {"uniprot", "uniprot_id", "entry"}:
                values.add(normalize_id(value))
    return values


def load_sources(args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, set[str]], dict[str, int]]:
    mcshane = pd.read_excel(args.mcshane, sheet_name="FinalTable_RPE_")
    mcshane = mcshane[["Protein IDs (Uniprot)", "Degradation profile"]].rename(
        columns={"Protein IDs (Uniprot)": "UniProt_ID", "Degradation profile": "Degradation_profile"}
    )
    mcshane["UniProt_ID"] = mcshane["UniProt_ID"].map(normalize_id)

    plddt = pd.read_csv(args.plddt)[["gene", "<pLDDT>", "std"]].rename(
        columns={"gene": "UniProt_ID", "<pLDDT>": "Average_pLDDT", "std": "pLDDT_SD"}
    )
    plddt["UniProt_ID"] = plddt["UniProt_ID"].map(normalize_id)
    plddt["Average_pLDDT"] = pd.to_numeric(plddt["Average_pLDDT"], errors="coerce")
    plddt["pLDDT_SD"] = pd.to_numeric(plddt["pLDDT_SD"], errors="coerce")
    if plddt["UniProt_ID"].duplicated().any():
        raise ValueError("Human_Avg_pLDDTs.csv contains duplicate UniProt IDs")

    entanglement = pd.read_csv(args.entanglement, sep="|", low_memory=False)[["gene", "CCBond"]]
    entanglement["gene"] = entanglement["gene"].map(normalize_id)
    entanglement = (
        entanglement.groupby("gene", as_index=False)
        .agg({"CCBond": merge_ccbond})
        .rename(columns={"gene": "UniProt_ID", "CCBond": "NCLE"})
    )

    knots = pd.read_table(args.knots, low_memory=False)
    knot_ids = set(
        knots.loc[knots["Organism"].eq("Human"), "Uniprot"]
        .astype(str)
        .str.strip()
        .str.split("-")
        .str[0]
    )
    proteome_lengths = parse_fasta(args.proteome_fasta)
    membrane_ids = set(parse_fasta(args.membrane_fasta))
    secreted_ids = parse_id_source(args.secreted_source)

    monomer_ids: set[str] = set()
    for path in args.monomer_fasta:
        monomer_ids.update(parse_fasta(path))

    master = mcshane.merge(plddt, on="UniProt_ID", how="inner", validate="one_to_one")
    master = master.loc[master["Average_pLDDT"] >= 70].copy()
    master = master.merge(entanglement, on="UniProt_ID", how="left", validate="one_to_one")
    master["Knotted"] = master["UniProt_ID"].isin(knot_ids).astype(int)
    master["FASTA_available"] = master["UniProt_ID"].isin(proteome_lengths).astype(int)
    master["Protein_length"] = master["UniProt_ID"].map(proteome_lengths)
    master["Membrane"] = master["UniProt_ID"].isin(membrane_ids).astype(int)
    master["Secreted"] = master["UniProt_ID"].isin(secreted_ids).astype(int)
    master["Known_monomer_annotation"] = np.where(
        master["UniProt_ID"].isin(monomer_ids), "Known_monomer", "Unclassified"
    )
    master["Entanglement_status"] = master["NCLE"].map(
        {1.0: "Native_entangled", 0.0: "Not_entangled"}
    ).fillna("Excluded_covalent_lasso")

    id_sets = {
        "knots": knot_ids,
        "membrane": membrane_ids,
        "secreted": secreted_ids,
        "known_monomer": monomer_ids,
    }
    return master, id_sets, proteome_lengths


def build_cohort(master: pd.DataFrame, profiles: list[str], name: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    flow: list[dict[str, object]] = []
    data = master.loc[master["Degradation_profile"].isin(profiles)].copy()

    def record(step: str) -> None:
        counts = data["Degradation_profile"].value_counts()
        flow.append(
            {
                "Analysis": name,
                "Step": step,
                "Total": int(len(data)),
                "ED": int(counts.get("ED", 0)),
                "NED": int(counts.get("NED", 0)),
                "UN": int(counts.get("UN", 0)),
            }
        )

    record("pLDDT >=70 and selected degradation categories")
    data = data.loc[data["NCLE"].notna()].copy()
    record("Usable binary NCLE status")
    data = data.loc[data["Knotted"].eq(0)].copy()
    record("Exclude knotted proteins")
    data = data.loc[data["FASTA_available"].eq(1)].copy()
    record("Require FASTA sequence")
    data = data.loc[data["Membrane"].eq(0)].copy()
    record("Exclude membrane proteins")
    data = data.loc[data["Secreted"].eq(0)].copy()
    record("Exclude secreted proteins")

    data["NED_binary"] = data["Degradation_profile"].eq("NED").astype(int)
    data["z_Protein_length"] = (
        data["Protein_length"] - data["Protein_length"].mean()
    ) / data["Protein_length"].std(ddof=1)
    data["Analysis_group"] = (
        data["Degradation_profile"]
        + "_"
        + data["NCLE"].map({1.0: "Entangled", 0.0: "NotEntangled"})
    )
    data = data.sort_values("UniProt_ID").reset_index(drop=True)
    data.insert(0, "Index", np.arange(1, len(data) + 1))
    return data, flow


def crude_analysis(data: pd.DataFrame, analysis: str) -> tuple[dict[str, object], pd.DataFrame]:
    a = int(((data["NCLE"] == 1) & (data["NED_binary"] == 1)).sum())
    b = int(((data["NCLE"] == 1) & (data["NED_binary"] == 0)).sum())
    c = int(((data["NCLE"] == 0) & (data["NED_binary"] == 1)).sum())
    d = int(((data["NCLE"] == 0) & (data["NED_binary"] == 0)).sum())
    table = np.array([[a, b], [c, d]])
    odds_ratio, fisher_p = stats.fisher_exact(table, alternative="two-sided")
    log_se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lower, upper = np.exp(np.log(odds_ratio) + np.array([-1.96, 1.96]) * log_se)
    result = {
        "Analysis": analysis,
        "Method": "Crude Fisher exact",
        "Term": "NCLE",
        "N": int(len(data)),
        "Coefficient": float(np.log(odds_ratio)),
        "Odds_ratio": float(odds_ratio),
        "CI_lower": float(lower),
        "CI_upper": float(upper),
        "P_value": float(fisher_p),
    }
    counts = pd.DataFrame(
        [
            {"Analysis": analysis, "NCLE_status": "Entangled", "NED": a, "Non_NED": b, "Total": a + b},
            {"Analysis": analysis, "NCLE_status": "Not_entangled", "NED": c, "Non_NED": d, "Total": c + d},
        ]
    )
    counts["NED_frequency"] = counts["NED"] / counts["Total"]
    return result, counts


def logistic_analysis(data: pd.DataFrame, analysis: str, method: str, terms: list[str]) -> list[dict[str, object]]:
    clean = data[["NED_binary", *terms]].dropna().copy()
    term_names = ["Intercept", *terms]
    x = np.column_stack([np.ones(len(clean)), clean[terms].to_numpy(dtype=float)])
    y = clean["NED_binary"].to_numpy(dtype=float)
    beta = np.zeros(x.shape[1], dtype=float)
    for _ in range(100):
        eta = np.clip(x @ beta, -35, 35)
        probability = 1.0 / (1.0 + np.exp(-eta))
        weights = np.clip(probability * (1.0 - probability), 1e-12, None)
        information = x.T @ (weights[:, None] * x)
        score = x.T @ (y - probability)
        step = np.linalg.solve(information, score)
        beta_next = beta + step
        if np.max(np.abs(step)) < 1e-10:
            beta = beta_next
            break
        beta = beta_next
    else:
        raise RuntimeError(f"Logistic regression did not converge for {analysis}: {method}")

    eta = np.clip(x @ beta, -35, 35)
    probability = 1.0 / (1.0 + np.exp(-eta))
    weights = np.clip(probability * (1.0 - probability), 1e-12, None)
    covariance = np.linalg.inv(x.T @ (weights[:, None] * x))
    standard_errors = np.sqrt(np.diag(covariance))
    coefficients = dict(zip(term_names, beta))
    errors = dict(zip(term_names, standard_errors))
    rows: list[dict[str, object]] = []
    for term in terms:
        coefficient = float(coefficients[term])
        standard_error = float(errors[term])
        lower_coefficient = coefficient - 1.96 * standard_error
        upper_coefficient = coefficient + 1.96 * standard_error
        rows.append(
            {
                "Analysis": analysis,
                "Method": method,
                "Term": term,
                "N": int(len(clean)),
                "Coefficient": coefficient,
                "Odds_ratio": float(np.exp(coefficient)),
                "CI_lower": float(np.exp(lower_coefficient)),
                "CI_upper": float(np.exp(upper_coefficient)),
                "P_value": float(2 * stats.norm.sf(abs(coefficient / standard_error))),
            }
        )
    return rows


def add_complex_annotation(data: pd.DataFrame, path: Path) -> pd.DataFrame:
    annotation = pd.read_csv(path)
    expected = {"UniProt_ID", "Complex_member"}
    if not expected.issubset(annotation.columns):
        raise ValueError(f"Complex annotation must contain columns: {sorted(expected)}")
    annotation = annotation[["UniProt_ID", "Complex_member"]].copy()
    annotation["UniProt_ID"] = annotation["UniProt_ID"].map(normalize_id)
    annotation["Complex_member"] = pd.to_numeric(annotation["Complex_member"], errors="coerce")
    if annotation["UniProt_ID"].duplicated().any():
        raise ValueError("Complex annotation contains duplicate UniProt IDs")
    if not set(annotation["Complex_member"].dropna().unique()).issubset({0, 1}):
        raise ValueError("Complex_member must be 0, 1, or missing")
    return data.merge(annotation, on="UniProt_ID", how="left", validate="one_to_one")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    master, id_sets, proteome_lengths = load_sources(args)

    primary, flow_primary = build_cohort(master, ["ED", "NED"], "Primary: NED vs ED")
    sensitivity, flow_sensitivity = build_cohort(
        master, ["ED", "NED", "UN"], "Sensitivity: NED vs ED + UN"
    )

    results: list[dict[str, object]] = []
    count_tables: list[pd.DataFrame] = []
    for label, data in [
        ("Primary: NED vs ED", primary),
        ("Sensitivity: NED vs ED + UN", sensitivity),
    ]:
        crude, counts = crude_analysis(data, label)
        results.append(crude)
        count_tables.append(counts)
        results.extend(
            logistic_analysis(
                data,
                label,
                "Model 1: NED ~ NCLE + z(length)",
                ["NCLE", "z_Protein_length"],
            )
        )

    complex_status = "Not run: accurate binary complex-membership annotation was not supplied."
    if args.complex_annotation:
        complex_status = "Run using the supplied Complex_member annotation."
        updated: list[pd.DataFrame] = []
        for label, data in [
            ("Primary: NED vs ED", primary),
            ("Sensitivity: NED vs ED + UN", sensitivity),
        ]:
            annotated = add_complex_annotation(data, args.complex_annotation)
            results.extend(
                logistic_analysis(
                    annotated,
                    label,
                    "Model 2: NED ~ NCLE + complex + z(length)",
                    ["NCLE", "Complex_member", "z_Protein_length"],
                )
            )
            updated.append(annotated)
        primary, sensitivity = updated

    results_df = pd.DataFrame(results)
    counts_df = pd.concat(count_tables, ignore_index=True)
    flow_df = pd.DataFrame(flow_primary + flow_sensitivity)

    ordered_columns = [
        "Index",
        "UniProt_ID",
        "Degradation_profile",
        "NED_binary",
        "NCLE",
        "Entanglement_status",
        "Average_pLDDT",
        "pLDDT_SD",
        "Protein_length",
        "z_Protein_length",
        "Knotted",
        "FASTA_available",
        "Membrane",
        "Secreted",
        "Known_monomer_annotation",
        "Analysis_group",
    ]
    if "Complex_member" in primary.columns:
        ordered_columns.append("Complex_member")
    primary[ordered_columns].to_csv(args.output_dir / "primary_cohort.csv", index=False)
    sensitivity[ordered_columns].to_csv(args.output_dir / "sensitivity_cohort.csv", index=False)
    results_df.to_csv(args.output_dir / "association_results.csv", index=False)
    counts_df.to_csv(args.output_dir / "contingency_counts.csv", index=False)
    flow_df.to_csv(args.output_dir / "cohort_flow.csv", index=False)

    metadata = {
        "primary_n": int(len(primary)),
        "sensitivity_n": int(len(sensitivity)),
        "primary_group_counts": primary["Analysis_group"].value_counts().sort_index().to_dict(),
        "sensitivity_profile_ncle_counts": {
            f"{profile}|NCLE={int(ncle)}": int(count)
            for (profile, ncle), count in sensitivity.groupby(
                ["Degradation_profile", "NCLE"]
            ).size().items()
        },
        "secreted_source_ids": len(id_sets["secreted"]),
        "known_monomer_ids": len(id_sets["known_monomer"]),
        "complex_model_status": complex_status,
        "notes": [
            "MobiDB disorder was not used to select either statistical cohort.",
            "NCLE=1 denotes native non-covalent entanglement; NCLE=0 denotes no detected entanglement.",
            "Proteins with covalent lasso records are excluded from binary NCLE analyses.",
            "Known_monomer_annotation=Unclassified must not be interpreted as complex membership.",
        ],
    }
    (args.output_dir / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2, default=str))

    print(json.dumps(metadata, indent=2, default=str))
    print("\nAssociation results:")
    print(results_df.to_string(index=False))


if __name__ == "__main__":
    main()
