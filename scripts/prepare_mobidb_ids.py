#!/usr/bin/env python3
"""Export query IDs using Jiang's historical average-pLDDT file.

This reconstructed utility prepares inputs; it does not query MobiDB.
The inclusive cutoff reproduces the 13,690-protein broad query list.
"""
import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plddt", type=Path,
                        default=ROOT / "raw_sources/Human_Avg_pLDDTs.csv")
    parser.add_argument("--out", type=Path, default=ROOT / "results/query_inputs")
    parser.add_argument("--cutoff", type=float, default=70)
    args = parser.parse_args()
    data = pd.read_csv(args.plddt)
    if data["gene"].isna().any() or data["gene"].duplicated().any():
        raise ValueError("Missing or repeated pLDDT accession")
    scores = pd.to_numeric(data["<pLDDT>"], errors="raise")
    selected = data.loc[scores.ge(args.cutoff)].copy()
    ids = selected["gene"].astype(str).str.strip()
    if ids.eq("").any() or ids.duplicated().any():
        raise ValueError("Invalid normalized query ID")
    args.out.mkdir(parents=True, exist_ok=True)
    cutoff_tag = f"{args.cutoff:g}".replace(".", "p")
    ids.to_csv(args.out / f"uniprots_plddt_ge{cutoff_tag}.txt", index=False, header=False)
    ids.head(5).to_csv(args.out / "uniprots_test5.txt", index=False, header=False)
    selected.to_csv(args.out / "selected_plddt_source.csv", index=False)
    print(f"Average pLDDT >= {args.cutoff:g}: {len(ids):,} unique proteins")
    print(f"Query inputs: {args.out.resolve()}")


if __name__ == "__main__":
    main()
