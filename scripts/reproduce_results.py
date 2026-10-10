#!/usr/bin/env python3
"""Recompute human results from frozen protein rows, not stored ORs.

This is a new reproduction runner. statistical_analysis.py is the recovered
original analysis code. This runner never refreshes online annotations.
"""
import argparse
import hashlib
import importlib.metadata
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import fisher_exact, norm

from statistical_analysis import crude_analysis, logistic_analysis

ROOT = Path(__file__).resolve().parents[1]
LABELS = {
    "Primary": "Primary: NED vs ED",
    "Sensitivity": "Sensitivity: NED vs ED + UN",
}
M1 = "Model 1: NED ~ NCLE + z(length)"
M2 = "Model 2: NED ~ NCLE + complex + z(length)"
METRICS = ["N", "Coefficient", "Odds_ratio", "CI_lower", "CI_upper", "P_value"]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_cohort(data, primary):
    if data["UniProt_ID"].isna().any() or data["UniProt_ID"].duplicated().any():
        raise ValueError("Missing or repeated UniProt ID")
    allowed = {"ED", "NED"} if primary else {"ED", "NED", "UN"}
    if not set(data["Degradation_profile"]).issubset(allowed):
        raise ValueError("Unexpected degradation category")
    for col in ["NCLE", "NED_binary"]:
        if data[col].isna().any() or not set(data[col]).issubset({0, 1}):
            raise ValueError(f"{col} must be complete and binary")
    expected_y = data["Degradation_profile"].eq("NED").astype(int)
    if not data["NED_binary"].eq(expected_y).all():
        raise ValueError("Outcome does not agree with degradation profile")
    length = pd.to_numeric(data["Protein_length"], errors="raise")
    if length.isna().any() or not length.gt(0).all() or length.std(ddof=1) <= 0:
        raise ValueError("Invalid protein lengths")
    if not data["Average_pLDDT"].ge(70).all():
        raise ValueError("A protein fails the inclusive pLDDT criterion")
    for col in ["Knotted", "Membrane", "Secreted"]:
        if not data[col].eq(0).all():
            raise ValueError(f"Excluded {col} protein in cohort")
    if not data["FASTA_available"].eq(1).all():
        raise ValueError("A FASTA sequence is missing")
    if not set(data["Complex_member"].dropna()).issubset({0, 1}):
        raise ValueError("Invalid complex-membership value")
    z = (length - length.mean()) / length.std(ddof=1)
    if not np.allclose(z, data["z_Protein_length"], atol=1e-12, rtol=1e-12):
        raise ValueError("Stored standardized lengths are inconsistent")
    # Standardize once within each full cohort, before Model 2 exclusions.
    data = data.copy()
    data["z_Protein_length"] = z
    return data


def independent_fit(data, label):
    """Independent MLE fit using a different optimizer and analytic Hessian."""
    terms = ["NCLE", "Complex_member", "z_Protein_length"]
    clean = data.dropna(subset=["NED_binary", *terms])
    x = np.column_stack([np.ones(len(clean)), clean[terms].to_numpy(float)])
    y = clean["NED_binary"].to_numpy(float)
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError("Model 2 design matrix is not full rank")

    def objective(beta):
        eta = x @ beta
        return np.logaddexp(0, eta).sum() - y @ eta

    def gradient(beta):
        return x.T @ (expit(x @ beta) - y)

    def hessian(beta):
        p = expit(x @ beta)
        return x.T @ ((p * (1 - p))[:, None] * x)

    fit = minimize(objective, np.zeros(x.shape[1]), jac=gradient,
                   hess=hessian, method="trust-exact",
                   options={"gtol": 1e-9, "maxiter": 500})
    if np.max(np.abs(gradient(fit.x))) >= 1e-6:
        raise RuntimeError(f"Independent fit did not converge: {fit.message}")
    se = np.sqrt(np.diag(np.linalg.inv(hessian(fit.x))))
    rows = []
    for i, term in enumerate(terms, 1):
        b = fit.x[i]
        rows.append(dict(Analysis=label, Method=M2, Term=term, N=len(clean),
                         Coefficient=b, Odds_ratio=np.exp(b),
                         CI_lower=np.exp(b - 1.96 * se[i]),
                         CI_upper=np.exp(b + 1.96 * se[i]),
                         P_value=2 * norm.sf(abs(b / se[i]))))
    return rows


def compare_results(actual, reference, description):
    keys = ["Analysis", "Method", "Term"]
    if actual.duplicated(keys).any() or reference.duplicated(keys).any():
        raise ValueError("Repeated result key")
    joined = actual.merge(reference, on=keys, how="outer",
                          suffixes=("_actual", "_reference"),
                          indicator=True, validate="one_to_one")
    rows = []
    for _, row in joined.iterrows():
        for metric in METRICS:
            a, b = row[f"{metric}_actual"], row[f"{metric}_reference"]
            match = row["_merge"] == "both" and np.isfinite(a) and np.isfinite(b)
            match = match and (a == b if metric == "N" else
                               np.isclose(a, b, atol=1e-8, rtol=1e-8))
            rows.append({**{k: row[k] for k in keys}, "Comparison": description,
                         "Metric": metric, "Actual": a, "Reference": b,
                         "Absolute_difference": abs(a - b), "Pass": bool(match)})
    return pd.DataFrame(rows)


def binary_diagnostic(data, exposed, comparison, label):
    y = data["NED_binary"].eq(1)
    a, b = int((exposed & y).sum()), int((exposed & ~y).sum())
    c, d = int((~exposed & y).sum()), int((~exposed & ~y).sum())
    odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return dict(Analysis=label, Comparison=comparison, Exposed_NED=a,
                Exposed_non_NED=b, Reference_NED=c, Reference_non_NED=d,
                Odds_ratio=float(odds), P_value=float(p))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path,
                        default=ROOT / "data/reference_human_analysis.xlsx")
    parser.add_argument("--out", type=Path, default=ROOT / "results/frozen_replay")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    cohorts = {}
    for name in LABELS:
        d = pd.read_excel(args.reference, sheet_name=f"{name} Cohort")
        cohorts[name] = validate_cohort(d, name == "Primary")
    primary_ids = set(cohorts["Primary"]["UniProt_ID"])
    sensitivity_edned = cohorts["Sensitivity"].query("Degradation_profile != 'UN'")
    if primary_ids != set(sensitivity_edned["UniProt_ID"]):
        raise ValueError("Primary and sensitivity ED/NED memberships differ")

    results, counts, restricted, exclusions, diagnostics = [], [], [], [], []
    independent, status_counts = [], []
    for name, data in cohorts.items():
        label = LABELS[name]
        crude, table = crude_analysis(data, label)
        results.append(crude)
        counts.append(table)
        results.extend(logistic_analysis(data, label, M1, ["NCLE", "z_Protein_length"]))
        results.extend(logistic_analysis(data, label, M2,
                                        ["NCLE", "Complex_member", "z_Protein_length"]))
        complete = data.dropna(subset=["Complex_member"])
        restricted.extend(logistic_analysis(complete, label,
                          "Restricted Model 1: NED ~ NCLE + z(length)",
                          ["NCLE", "z_Protein_length"]))
        provisional = data["Manual_review_decision"].eq("Complex_or_oligomer")
        exclusions.extend(logistic_analysis(data.loc[~provisional], label,
                          "Model 2 excluding provisional positive calls",
                          ["NCLE", "Complex_member", "z_Protein_length"]))
        diagnostics.append(binary_diagnostic(complete, complete["Complex_member"].eq(1),
                           "Complex member versus explicit monomer", label))
        diagnostics.append(binary_diagnostic(data, data["Complex_member"].notna(),
                           "Usable versus missing complex annotation", label))
        for status, group in data.groupby("Oligomer_status", dropna=False):
            status_counts.append(dict(Analysis=label, Status=status, N=len(group),
                                      NED=int(group["NED_binary"].sum())))
        independent.extend(independent_fit(data, label))
        data.to_csv(args.out / f"{name.lower()}_cohort.csv", index=False)

    actual = pd.DataFrame(results)
    expected = pd.read_excel(args.reference, sheet_name="Association Results")
    comparison = compare_results(actual, expected, "Reported workbook results")
    independent_comparison = compare_results(pd.DataFrame(independent),
                              actual.loc[actual["Method"].eq(M2)],
                              "Independent optimizer versus original solver")
    outputs = {
        "association_results.csv": actual,
        "contingency_counts.csv": pd.concat(counts, ignore_index=True),
        "restricted_model1.csv": pd.DataFrame(restricted),
        "provisional_exclusion_results.csv": pd.DataFrame(exclusions),
        "model2_fisher_diagnostics.csv": pd.DataFrame(diagnostics),
        "annotation_status_counts.csv": pd.DataFrame(status_counts),
        "reference_comparison.csv": comparison,
        "independent_optimizer_comparison.csv": independent_comparison,
    }
    for name, data in outputs.items():
        data.to_csv(args.out / name, index=False)
    passed = comparison["Pass"].all() and independent_comparison["Pass"].all()
    manifest = {
        "run_utc": datetime.now(timezone.utc).isoformat(), "python": sys.version,
        "reference_sha256": sha256(args.reference),
        "original_analysis_sha256": sha256(Path(__file__).with_name("statistical_analysis.py")),
        "runner_sha256": sha256(__file__), "numerical_comparison_passed": bool(passed),
        "primary_n": len(cohorts["Primary"]),
        "sensitivity_n": len(cohorts["Sensitivity"]),
        "versions": {p: importlib.metadata.version(p)
                     for p in ["numpy", "pandas", "scipy", "openpyxl"]},
        "scope": "Frozen-cohort numerical replay; no online annotation refresh",
        "confidence_interval_multiplier": 1.96,
    }
    (args.out / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(actual.loc[actual["Term"].eq("NCLE")].to_string(index=False))
    print(f"\nNumerical reproduction: {'PASS' if passed else 'FAIL'}")
    print(f"Results: {args.out.resolve()}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
