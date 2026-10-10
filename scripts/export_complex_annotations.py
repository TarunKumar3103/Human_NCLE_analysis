#!/usr/bin/env python3
"""Export evidence and replay recorded review decisions from the final workbook.

This does not recreate the unavailable original text-classification code or
independently validate its biological decisions. It reconstructs final calls
from frozen automated calls plus the recorded first-pass review decisions.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path,
                        default=ROOT / "data/reference_human_analysis.xlsx")
    parser.add_argument("--out", type=Path, default=ROOT / "data")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    evidence = pd.read_excel(args.reference, sheet_name="UniProt Evidence")
    decisions = pd.read_excel(args.reference, sheet_name="Manual Review")
    if evidence["UniProt_ID"].isna().any() or evidence["UniProt_ID"].duplicated().any():
        raise ValueError("Evidence requires one row per protein")
    if decisions["UniProt_ID"].isna().any() or decisions["UniProt_ID"].duplicated().any():
        raise ValueError("Review decisions contain missing or repeated IDs")
    ambiguous_ids = set(evidence.loc[evidence["Automated_Oligomer_status"].eq("Ambiguous"),
                                     "UniProt_ID"])
    if set(decisions["UniProt_ID"]) != ambiguous_ids:
        raise ValueError("Review worksheet does not cover exactly the ambiguous records")
    decision_map = decisions.set_index("UniProt_ID")["Manual_review_decision"]
    if not set(decision_map.dropna()).issubset({"Complex_or_oligomer", "Unresolved"}):
        raise ValueError("Unexpected recorded review decision")
    recorded = evidence["UniProt_ID"].map(decision_map)
    if not recorded.fillna("").eq(evidence["Manual_review_decision"].fillna("")).all():
        raise ValueError("Evidence and review worksheets disagree")
    calls = evidence["Automated_Complex_member"].copy()
    # A workflow-completion flag never supplies a biological classification.
    positive = recorded.eq("Complex_or_oligomer")
    unresolved = recorded.eq("Unresolved")
    calls.loc[positive] = 1
    calls.loc[unresolved] = np.nan
    if not np.allclose(calls, evidence["Complex_member"], equal_nan=True):
        raise ValueError("Recorded review decisions do not reconstruct final calls")
    annotation = pd.DataFrame({"UniProt_ID": evidence["UniProt_ID"],
                               "Complex_member": calls})
    annotation.to_csv(args.out / "complex_annotations.csv", index=False)
    evidence.to_csv(args.out / "uniprot_evidence_snapshot.csv", index=False)
    decisions.to_csv(args.out / "complex_review_decisions.csv", index=False)
    print(f"Evidence proteins: {len(evidence)}")
    print(f"Initially ambiguous: {len(decisions)}")
    print(f"Provisional positive calls: {int(positive.sum())}")
    print(f"Reviewed unresolved: {int(unresolved.sum())}")
    print("Recorded-decision reconstruction: PASS")


if __name__ == "__main__":
    main()
