#!/usr/bin/env python3
"""Audit the frozen processed MobiDB export and its statistical-cohort coverage.

The export is processed output, not raw API JSON. Coordinate clipping is
reported explicitly; it does not establish sequence-version compatibility.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def intervals(value):
    if pd.isna(value) or str(value).strip() in {"", "-"}:
        return []
    parsed = []
    for token in str(value).split(";"):
        start, end = map(int, token.strip().split("-"))
        if start > end:
            raise ValueError(f"Reversed interval: {token}")
        parsed.append((start, end))
    parsed.sort()
    merged = []
    for start, end in parsed:
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master", type=Path, default=ROOT / "data/processed_MobiDB_master.csv")
    parser.add_argument("--reference", type=Path, default=ROOT / "data/reference_human_analysis.xlsx")
    parser.add_argument("--out", type=Path, default=ROOT / "results/disorder_audit")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(args.master)
    if d["UniProt_ID"].isna().any() or d["UniProt_ID"].duplicated().any():
        raise ValueError("Master identifiers must be present and unique")
    fraction = pd.to_numeric(d["Fraction_IDP_IDR"], errors="raise")
    length = pd.to_numeric(d["MobiDB_length"], errors="raise")
    count = pd.to_numeric(d["Disorder_residues"], errors="raise")
    available = fraction.notna()
    if not fraction.loc[available].between(0, 1).all():
        raise ValueError("Out-of-range fraction")
    if not length.loc[available].gt(0).all() or not count.loc[available].ge(0).all():
        raise ValueError("Invalid denominator or disorder count")
    if not np.allclose(fraction.loc[available], (count / length).loc[available], atol=1e-8, rtol=0):
        raise ValueError("Fraction disagrees with count divided by length")
    ok = d["Disorder_query_status"].eq("ok")
    if not ok.eq(available).all():
        raise ValueError("Status and fraction availability disagree")
    issues = []
    for _, row in d.loc[available].iterrows():
        regions = intervals(row["Disorder_regions"])
        full = sum(end - start + 1 for start, end in regions)
        size = int(row["MobiDB_length"])
        clipped = sum(max(0, min(size, end) - max(1, start) + 1) for start, end in regions)
        if clipped != int(row["Disorder_residues"]):
            raise ValueError(f"Clipped count does not reproduce stored count: {row['UniProt_ID']}")
        if full != clipped:
            issues.append(dict(UniProt_ID=row["UniProt_ID"], Full_interval_count=full,
                               Clipped_count=clipped, MobiDB_length=size,
                               Disorder_regions=row["Disorder_regions"]))
    pd.DataFrame(issues).to_csv(args.out / "interval_coordinate_warnings.csv", index=False)
    rows = []
    for name in ["Primary", "Sensitivity"]:
        cohort = pd.read_excel(args.reference, sheet_name=f"{name} Cohort")
        merged = cohort.merge(d[["UniProt_ID", "Fraction_IDP_IDR"]],
                              how="left", on="UniProt_ID", validate="one_to_one")
        merged["Fraction_available"] = merged["Fraction_IDP_IDR"].notna()
        for group, sub in merged.groupby("Analysis_group"):
            n, found = len(sub), int(sub["Fraction_available"].sum())
            rows.append(dict(Analysis=name, Group=group, Total=n, Fraction_available=found,
                             Fraction_missing=n-found, Missing_percent=100*(n-found)/n))
        # Coverage output has no disorder cutoff and is not the statistical sample.
        merged.to_csv(args.out / f"{name.lower()}_cohort_with_disorder.csv", index=False)
    pd.DataFrame(rows).to_csv(args.out / "disorder_coverage.csv", index=False)
    summary = dict(master_n=len(d), status_counts=d["Disorder_query_status"].value_counts().to_dict(),
                   fraction_available=int(available.sum()), interval_coordinate_warnings=len(issues),
                   source_kind="Processed export; raw API responses are not included")
    (args.out / "disorder_audit_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
