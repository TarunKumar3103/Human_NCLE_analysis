#!/usr/bin/env python3
"""Replay raw cohort filters and compare them with the final reference.

Without --secreted-source, validation stops before secreted exclusion. The
final eligible IDs are never used to manufacture a secreted annotation list.
"""
import argparse
import hashlib
import json
import re
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from statistical_analysis import load_sources, build_cohort

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, default=ROOT / "raw_sources")
    parser.add_argument("--reference", type=Path,
                        default=ROOT / "data/reference_human_analysis.xlsx")
    parser.add_argument("--secreted-source", type=Path)
    parser.add_argument("--out", type=Path, default=ROOT / "results/raw_source_audit")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    sources = {
        "mcshane": args.sources / "McShane_human_mmc4.xlsx",
        "entanglement": args.sources / "Human_AF_combined_20250614.csv",
        "plddt": args.sources / "Human_Avg_pLDDTs.csv",
        "knots": args.sources / "AF4_knots.dat",
        "proteome_fasta": args.sources / "UP000005640_9606.fasta",
        "membrane_fasta": args.sources / "uniprotkb_human_keywords_membrane_2025_03_19.fasta",
    }
    missing = [str(p) for p in sources.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing raw source files: " + ", ".join(missing))
    for key in ["proteome_fasta", "membrane_fasta"]:
        ids = []
        with sources[key].open() as handle:
            for line in handle:
                if line.startswith(">"):
                    pieces = line[1:].strip().split("|")
                    ids.append(pieces[1] if len(pieces) >= 3 else pieces[0].split()[0])
        if len(ids) != len(set(ids)):
            raise ValueError(f"Duplicate accessions in {key}; do not concatenate sequences")
    full = pd.read_excel(sources["mcshane"], sheet_name="FinalTable_RPE_")
    if full["Protein IDs (Uniprot)"].isna().any() or full["Protein IDs (Uniprot)"].duplicated().any():
        raise ValueError("McShane identifiers are missing or duplicated")
    ent = pd.read_csv(sources["entanglement"], sep="|", usecols=["gene", "CCBond"])
    tokens = set(ent["CCBond"].dropna().astype(str).str.strip())
    if tokens - {"True", "False", "1", "0", "1.0", "0.0"}:
        raise ValueError(f"Unknown CCBond tokens: {tokens}")
    with tempfile.TemporaryDirectory() as tmp:
        empty = Path(tmp) / "empty_for_presecreted_validation.json"
        empty.write_text("[]")
        ns = SimpleNamespace(**sources, secreted_source=args.secreted_source or empty,
                             monomer_fasta=[])
        master, sets, lengths = load_sources(ns)
    layout = pd.read_excel(args.reference, sheet_name="Cohort Flow", header=None)
    comparisons, membership = [], []
    for name, profiles in [("Primary", ["ED", "NED"]),
                           ("Sensitivity", ["ED", "NED", "UN"])]:
        label = "Primary: NED vs ED" if name == "Primary" else "Sensitivity: NED vs ED + UN"
        data, flow = build_cohort(master, profiles, label)
        if args.secreted_source is None:
            flow = flow[:-1]
            data["Secreted"] = np.nan
        for entry in flow:
            ref = layout.loc[layout[0].eq(label) & layout[1].eq(entry["Step"])]
            if len(ref) != 1:
                raise ValueError("Reference flow step is missing or repeated")
            for i, metric in enumerate(["Total", "ED", "NED", "UN"], 2):
                expected = int(ref.iloc[0, i])
                comparisons.append(dict(Analysis=label, Step=entry["Step"], Metric=metric,
                                        Actual=entry[metric], Reference=expected,
                                        Pass=entry[metric] == expected))
        prefix = "rebuilt" if args.secreted_source else "pre_secreted"
        data.to_csv(args.out / f"{prefix}_{name.lower()}_cohort.csv", index=False)
        if args.secreted_source:
            ref = pd.read_excel(args.reference, sheet_name=f"{name} Cohort")
            joined = data.merge(ref, on="UniProt_ID", how="outer", indicator=True,
                                suffixes=("_rebuilt", "_reference"), validate="one_to_one")
            for _, row in joined.iterrows():
                same = row["_merge"] == "both"
                for field in ["NCLE", "NED_binary", "Protein_length", "Average_pLDDT"]:
                    same = same and np.isclose(row[field + "_rebuilt"], row[field + "_reference"],
                                               atol=1e-10, rtol=1e-10)
                membership.append(dict(Analysis=label, UniProt_ID=row["UniProt_ID"],
                                       Membership=row["_merge"], Pass=bool(same)))
    result = pd.DataFrame(comparisons)
    result.to_csv(args.out / "cohort_flow_comparison.csv", index=False)
    if membership:
        pd.DataFrame(membership).to_csv(args.out / "final_membership_comparison.csv", index=False)
    # Preserve the historical normalization, and expose its broader scope.
    knot = pd.read_table(sources["knots"])
    human = knot.loc[knot["Organism"].eq("Human"), "Uniprot"].astype(str).str.strip()
    broad = set(human.str.split("-").str[0])
    strict = {re.sub(r"-F\d+$", "", x) for x in human}
    qc_ids = set(pd.read_csv(sources["plddt"]).query("`<pLDDT>` >= 70")["gene"])
    meta = dict(mcshane_rows=len(full), degradation_counts=full["Degradation profile"].value_counts().to_dict(),
                entanglement_rows=len(ent), entanglement_proteins=ent["gene"].nunique(),
                historical_knotted_in_plddt_master=len(qc_ids & broad),
                strict_fragment_knotted_in_plddt_master=len(qc_ids & strict),
                broader_knot_normalization_extra_ids=sorted((qc_ids & broad) - strict),
                secreted_source_supplied=args.secreted_source is not None,
                final_cohort_source_validation="RUN" if args.secreted_source else "NOT RUN",
                source_sha256={k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in sources.items()})
    (args.out / "raw_source_manifest.json").write_text(json.dumps(meta, indent=2))
    print(result.to_string(index=False))
    print("Final secreted exclusion: " + ("checked" if args.secreted_source else
          "NOT independently replayed; supply the approved historical keyword export"))
    if not result["Pass"].all() or (membership and not all(x["Pass"] for x in membership)):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
