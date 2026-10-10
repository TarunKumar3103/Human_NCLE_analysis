import time
import json
import argparse

import requests
from bioservices import UniProt


# MobiDB keywords in priority order
DISORDER_KEYWORDS = [
    'curated-disorder-priority',
    'homology-disorder-priority',
    'derived-missing_residues-priority',
    'prediction-disorder-priority',
    'derived-mobile-th_90',
]

MOBIDB_API_RATE_LIMIT_SEC = 0.5  # seconds between API calls


def get_mobidb_disordered_data(uniprot: str) -> list:
    """
    Retrieve disordered regions for a protein from MobiDB.

    Queries the MobiDB API for disorder annotations and returns a list of
    interval tuples ordered by the priority defined in DISORDER_KEYWORDS.

    Args:
        uniprot: UniProt accession ID.

    Returns:
        List of interval tuples, e.g. [([1, 56],), ([57, 68],)].
        Returns an empty list if the ID is not found or has no disorder data.
    """
    url = f'https://mobidb.org/api/download?format=json&acc={uniprot}'

    try:
        res = requests.get(url, timeout=30)
        res.raise_for_status()
    except requests.RequestException as e:
        print(f"[WARNING] MobiDB request failed for {uniprot}: {e}")
        return []

    try:
        result = res.json()
    except json.JSONDecodeError:
        print(f"[WARNING] Could not parse MobiDB JSON response for {uniprot}")
        return []

    disordered_regions = []
    for key in DISORDER_KEYWORDS:
        if key in result:
            regions = result[key]['regions']
            disordered_regions.append(tuple(regions))

    return disordered_regions


def merge_data(data_sets: list) -> tuple:
    """
    Merge and compact intervals from multiple disorder data sets.

    Combines all intervals, sorts by start position, merges overlapping
    intervals, and coalesces adjacent intervals.

    Args:
        data_sets: List of interval tuples as returned by get_mobidb_disordered_data.

    Returns:
        Tuple of merged [start, end] intervals, or an empty tuple if no
        valid intervals are present.
    """
    if not data_sets or isinstance(data_sets, str):
        return ()

    all_intervals = []
    for intervals in data_sets:
        all_intervals.extend(intervals)

    if not all_intervals:
        return ()

    # Sort by start position
    all_intervals.sort(key=lambda x: x[0])

    # Merge overlapping intervals
    merged = [list(all_intervals[0])]
    for interval in all_intervals[1:]:
        if interval[0] <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], interval[1])
        else:
            merged.append(list(interval))

    # Coalesce adjacent intervals (end + 1 == next start)
    final = [merged[0]]
    for interval in merged[1:]:
        if final[-1][1] + 1 == interval[0]:
            final[-1][1] = interval[1]
        else:
            final.append(interval)

    return tuple(final)


def count_overlap(first_range: tuple, second_ranges: tuple) -> list:
    """
    Count overlapping residues between a query range and a set of disorder regions.

    Args:
        first_range: (start, end) tuple for the query range (e.g. full protein length).
        second_ranges: Tuple of [start, end] lists representing disorder regions.

    Returns:
        List of overlap counts per region. Returns an empty list if second_ranges is empty.
    """
    if not second_ranges:
        return []

    overlaps = []
    for region in second_ranges:
        start = max(first_range[0], region[0])
        end = min(first_range[1], region[1])
        overlaps.append(max(0, end - start + 1))

    return overlaps


def get_uniprot_length(uniprot: str, u: UniProt) -> int | None:
    """
    Retrieve the sequence length of a protein from UniProt.

    Args:
        uniprot: UniProt accession ID.
        u: Initialised bioservices UniProt instance.

    Returns:
        Sequence length as an integer, or None if the lookup fails.
    """
    res = u.search(uniprot, frmt="tsv", columns="length", limit=1, database='uniprotkb')
    lines = res.strip().split('\n') if res else []

    if len(lines) < 2:
        print(f"[WARNING] Could not find length for {uniprot} in UniProt")
        return None

    try:
        return int(lines[1])
    except ValueError:
        print(f"[WARNING] Unexpected length value for {uniprot}: {lines[1]!r}")
        return None


def get_disorder_summary(uniprot: str, u: UniProt) -> dict:
    """
    Compute the disorder summary for a single protein.

    Args:
        uniprot: UniProt accession ID.
        u: Initialised bioservices UniProt instance.

    Returns:
        Dict with keys:
            'uniprot'                : accession ID
            'length'                 : protein length (int) or None on failure
            'disorder_regions'       : tuple of merged [start, end] disorder intervals
            'num_disorder_residues'  : total number of disordered residues (int) or None
    """
    length = get_uniprot_length(uniprot, u)
    if length is None:
        return {
            'uniprot': uniprot,
            'length': None,
            'disorder_regions': (),
            'num_disorder_residues': None,
        }

    disorder_regions = merge_data(get_mobidb_disordered_data(uniprot))
    overlap_counts = count_overlap((1, length), disorder_regions)
    num_disorder_residues = sum(overlap_counts)

    return {
        'uniprot': uniprot,
        'length': length,
        'disorder_regions': disorder_regions,
        'num_disorder_residues': num_disorder_residues,
    }


def format_regions(disorder_regions: tuple) -> str:
    """Render merged disorder intervals as a compact string, e.g. '1-56;70-88'."""
    if not disorder_regions:
        return '-'
    return ';'.join(f"{start}-{end}" for start, end in disorder_regions)


def read_uniprot_list(path: str) -> list:
    """Read UniProt accessions from a file, one per line, skipping blank lines."""
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def process_uniprots(uniprots: list, u: UniProt, output_file: str | None = None) -> None:
    """
    Fetch and report the disorder summary for each UniProt accession.

    Writes a tab-separated table (header + one row per accession) to stdout,
    and to output_file as well when provided. Columns:
        uniprot  disorder_regions  num_disorder_residues  length

    Args:
        uniprots: List of UniProt accession IDs.
        u: Initialised bioservices UniProt instance.
        output_file: Optional path to also write the table to.
    """
    header = "uniprot\tdisorder_regions\tnum_disorder_residues\tlength"

    out = open(output_file, 'w') if output_file else None
    try:
        print(header)
        if out:
            out.write(header + '\n')
            out.flush()

        for i, uniprot in enumerate(uniprots):
            summary = get_disorder_summary(uniprot, u)
            line = (
                f"{summary['uniprot']}\t"
                f"{format_regions(summary['disorder_regions'])}\t"
                f"{summary['num_disorder_residues']}\t"
                f"{summary['length']}"
            )
            print(line)
            if out:
                out.write(line + '\n')
                out.flush()

            # Rate-limit only between MobiDB API calls, not after the last one.
            if i < len(uniprots) - 1:
                time.sleep(MOBIDB_API_RATE_LIMIT_SEC)
    finally:
        if out:
            out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch disorder statistics (merged disorder range, number of "
                    "disordered residues, protein length) from MobiDB/UniProt for "
                    "one or more UniProt accessions."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("-u", "--uniprot",
                        help="A single UniProt accession ID.")
    source.add_argument("-f", "--file",
                        help="Path to a file listing UniProt accession IDs, one per line.")
    parser.add_argument("-o", "--output", default=None,
                        help="Optional output file path (tab-separated). "
                             "Results are always printed to stdout.")
    args = parser.parse_args()

    if args.uniprot:
        uniprots = [args.uniprot.strip()]
    else:
        uniprots = read_uniprot_list(args.file)

    if not uniprots:
        parser.error("No UniProt accessions found in the provided input.")

    u = UniProt()
    process_uniprots(uniprots, u, output_file=args.output)
