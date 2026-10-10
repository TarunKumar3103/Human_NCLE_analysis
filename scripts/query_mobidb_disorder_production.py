#!/usr/bin/env python3
"""Query MobiDB/UniProt and calculate the per-protein disorder fraction.

This production version preserves the disorder-track union used in Quyen's
sample script, while adding retries, explicit failure statuses, incremental
output, and resume support. Failed requests remain missing rather than being
misclassified as zero disorder.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Iterable
from urllib.parse import quote

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


DISORDER_KEYWORDS = [
    "curated-disorder-priority",
    "homology-disorder-priority",
    "derived-missing_residues-priority",
    "prediction-disorder-priority",
    "derived-mobile-th_90",
]

FIELDNAMES = [
    "uniprot",
    "disorder_regions",
    "num_disorder_residues",
    "length",
    "fraction_idp_idr",
    "query_status",
    "annotation_tracks",
]


def build_session() -> requests.Session:
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        status=5,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "protein-decay-disorder-table/1.0",
            "Accept": "application/json",
        }
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def request_json(session: requests.Session, url: str, timeout: int = 45):
    try:
        response = session.get(url, timeout=timeout)
    except requests.RequestException as exc:
        return None, f"request_failed:{type(exc).__name__}"

    if response.status_code == 404:
        return None, "not_found"
    if not response.ok:
        return None, f"http_{response.status_code}"

    try:
        return response.json(), None
    except (json.JSONDecodeError, ValueError):
        return None, "invalid_json"


def normalize_mobidb_payload(payload):
    if isinstance(payload, dict):
        return payload
    if isinstance(payload, list) and len(payload) == 1 and isinstance(payload[0], dict):
        return payload[0]
    return None


def get_mobidb_regions(session: requests.Session, uniprot: str):
    url = f"https://mobidb.org/api/download?format=json&acc={quote(uniprot)}"
    payload, error = request_json(session, url)
    if error:
        return None, [], f"mobidb_{error}"

    result = normalize_mobidb_payload(payload)
    if result is None:
        return None, [], "mobidb_unexpected_payload"

    data_sets = []
    tracks = []
    for key in DISORDER_KEYWORDS:
        block = result.get(key)
        if not isinstance(block, dict):
            continue
        regions = block.get("regions")
        if not isinstance(regions, list):
            continue
        tracks.append(key)
        data_sets.append(regions)

    if not tracks:
        return None, [], "mobidb_no_selected_annotation"

    return data_sets, tracks, None


def get_uniprot_length(session: requests.Session, uniprot: str):
    url = f"https://rest.uniprot.org/uniprotkb/{quote(uniprot)}.json"
    payload, error = request_json(session, url)
    if error:
        return None, f"uniprot_{error}"
    try:
        length = int(payload["sequence"]["length"])
    except (KeyError, TypeError, ValueError):
        return None, "uniprot_missing_length"
    if length <= 0:
        return None, "uniprot_invalid_length"
    return length, None


def merge_intervals(data_sets: Iterable[Iterable]) -> tuple[tuple[int, int], ...]:
    intervals = []
    for data_set in data_sets:
        for region in data_set:
            if not isinstance(region, (list, tuple)) or len(region) < 2:
                continue
            try:
                start = int(region[0])
                end = int(region[1])
            except (TypeError, ValueError):
                continue
            if start > end:
                start, end = end, start
            intervals.append((start, end))

    if not intervals:
        return ()

    intervals.sort(key=lambda item: (item[0], item[1]))
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)


def count_disordered_residues(regions, length: int) -> int:
    total = 0
    for start, end in regions:
        clipped_start = max(1, start)
        clipped_end = min(length, end)
        if clipped_end >= clipped_start:
            total += clipped_end - clipped_start + 1
    return total


def format_regions(regions) -> str:
    return ";".join(f"{start}-{end}" for start, end in regions) if regions else "-"


def query_one(session: requests.Session, uniprot: str) -> dict[str, str]:
    data_sets, tracks, mobidb_error = get_mobidb_regions(session, uniprot)
    length, length_error = get_uniprot_length(session, uniprot)

    record = {field: "" for field in FIELDNAMES}
    record["uniprot"] = uniprot
    record["annotation_tracks"] = ";".join(tracks)

    errors = [error for error in (mobidb_error, length_error) if error]
    if errors:
        if length is not None:
            record["length"] = str(length)
        record["query_status"] = "|".join(errors)
        return record

    regions = merge_intervals(data_sets)
    disorder_residues = count_disordered_residues(regions, length)
    fraction = disorder_residues / length
    if not 0.0 <= fraction <= 1.0:
        record["query_status"] = "invalid_fraction"
        return record

    record.update(
        {
            "disorder_regions": format_regions(regions),
            "num_disorder_residues": str(disorder_residues),
            "length": str(length),
            "fraction_idp_idr": f"{fraction:.8f}",
            "query_status": "ok",
        }
    )
    return record


def read_ids(path: Path) -> list[str]:
    seen = set()
    ids = []
    with path.open() as handle:
        for line in handle:
            accession = line.strip()
            if accession and accession not in seen:
                seen.add(accession)
                ids.append(accession)
    return ids


def read_existing(path: Path) -> OrderedDict[str, dict[str, str]]:
    records = OrderedDict()
    if not path.exists() or path.stat().st_size == 0:
        return records
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != FIELDNAMES:
            raise ValueError(
                f"Existing output has unexpected columns: {reader.fieldnames}. "
                "Use a new output filename."
            )
        for row in reader:
            accession = row.get("uniprot", "").strip()
            if accession:
                records[accession] = row
    return records


def compact_existing(path: Path, records, input_ids: list[str]) -> None:
    if not records:
        return
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        written = set()
        for accession in input_ids:
            if accession in records:
                writer.writerow(records[accession])
                written.add(accession)
        for accession, record in records.items():
            if accession not in written:
                writer.writerow(record)
    os.replace(tmp_path, path)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Query disorder annotations and calculate fraction_IDP_IDR."
    )
    parser.add_argument("-f", "--file", required=True, help="One UniProt ID per line.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path.")
    parser.add_argument(
        "--rate-limit",
        type=float,
        default=0.5,
        help="Seconds to pause after each protein (default: 0.5).",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help="Print progress every N processed proteins (default: 100).",
    )
    args = parser.parse_args()

    input_path = Path(args.file)
    output_path = Path(args.output)
    uniprots = read_ids(input_path)
    if not uniprots:
        parser.error("No UniProt IDs found in the input file.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    existing = read_existing(output_path)
    compact_existing(output_path, existing, uniprots)
    completed = {
        accession
        for accession, record in existing.items()
        if record.get("query_status") == "ok"
    }
    pending = [accession for accession in uniprots if accession not in completed]

    output_exists = output_path.exists() and output_path.stat().st_size > 0
    session = build_session()
    processed = 0
    print(
        f"Input IDs: {len(uniprots):,}; already successful: {len(completed):,}; "
        f"to query: {len(pending):,}",
        flush=True,
    )

    with output_path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        if not output_exists:
            writer.writeheader()
            handle.flush()

        try:
            for accession in pending:
                record = query_one(session, accession)
                writer.writerow(record)
                handle.flush()
                processed += 1

                if record["query_status"] != "ok":
                    print(
                        f"[WARNING] {accession}: {record['query_status']}",
                        file=sys.stderr,
                        flush=True,
                    )
                if args.progress_every > 0 and processed % args.progress_every == 0:
                    print(
                        f"Processed this run: {processed:,}/{len(pending):,}",
                        flush=True,
                    )
                if args.rate_limit > 0:
                    time.sleep(args.rate_limit)
        except KeyboardInterrupt:
            print("\nStopped safely. Rerun the same command to resume.", file=sys.stderr)
            return 130

    final_records = read_existing(output_path)
    compact_existing(output_path, final_records, uniprots)
    final_records = read_existing(output_path)
    counts = Counter(record.get("query_status", "") for record in final_records.values())
    print(f"Rows in output: {len(final_records):,}")
    print("Status counts:")
    for status, count in sorted(counts.items()):
        print(f"  {status or '<blank>'}: {count:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
