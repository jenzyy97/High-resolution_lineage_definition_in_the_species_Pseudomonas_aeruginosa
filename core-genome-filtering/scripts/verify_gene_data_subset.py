#!/usr/bin/env python3
"""
verify_gene_data_subset.py

Validate a subset gene_data.csv produced by subset_gene_data.py against the
core p/a it was built from. Checks:

  1. COVERAGE  - every gene id referenced in the p/a cells has a row in the
                 subset (a truncated/aborted run shows up here as missing ids).
  2. SEQUENCES - no row has an empty / NA sequence in the sequence column.
  3. INTEGRITY - every data row has the expected number of fields (a run killed
                 mid-write usually leaves a final partial row).
  4. DUPLICATES - reports ids appearing more than once (informational).

Exit code is non-zero if any hard check (coverage / sequences / integrity)
fails, so it can gate a pipeline.

Example
-------
  python verify_gene_data_subset.py \
      --pa core_gene_presence_absence.csv \
      --subset core4079_gene_data.csv
"""

from __future__ import annotations

import argparse
import csv
import sys

import pandas as pd

META_COLS = ["Gene", "Non-unique Gene name", "Annotation"]
CELL_SPLIT = ";"


def genome_cols(pa: pd.DataFrame) -> list:
    return [c for c in pa.columns if c not in META_COLS]


def wanted_ids(pa: pd.DataFrame) -> set:
    gcols = genome_cols(pa)
    flat = pd.Series(pa[gcols].to_numpy().ravel(), dtype="object").dropna()
    flat = flat[flat.astype(str).str.strip() != ""]
    if flat.empty:
        return set()
    parts = flat.astype(str).str.split(CELL_SPLIT).explode().str.strip()
    return set(parts[parts != ""].unique())


def run(args: argparse.Namespace) -> None:
    print(f"Loading core p/a: {args.pa} ...", flush=True)
    pa = pd.read_csv(args.pa, dtype=str)
    want = wanted_ids(pa)
    print(f"  {len(want):,} unique gene ids expected from the p/a.")

    print(f"Scanning subset (streaming): {args.subset} ...", flush=True)
    seen: set = set()
    dup = 0
    n_rows = 0
    n_empty_seq = 0
    bad_rows = 0
    header = None
    id_idx = seq_idx = None
    ncols = None

    # stream with the csv module so a truncated final line is caught, not silently parsed
    with open(args.subset, newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            sys.exit("subset file is empty (no header).")
        ncols = len(header)
        if args.id_col not in header:
            sys.exit(f"id column '{args.id_col}' not in subset header: {header}")
        if args.seq_col not in header:
            sys.exit(f"seq column '{args.seq_col}' not in subset header: {header}")
        id_idx = header.index(args.id_col)
        seq_idx = header.index(args.seq_col)

        for row in reader:
            n_rows += 1
            if len(row) != ncols:
                bad_rows += 1
                continue  # malformed / truncated row - don't trust its fields
            gid = row[id_idx].strip()
            seq = row[seq_idx].strip()
            if gid in seen:
                dup += 1
            else:
                seen.add(gid)
            if seq == "":
                n_empty_seq += 1
            if args.progress and n_rows % 5_000_000 == 0:
                print(f"  ... scanned {n_rows:,} rows", flush=True)

    missing = want - seen
    extra = seen - want  # ids in subset not referenced by the p/a (shouldn't happen)

    print("\n================ VERIFY REPORT ================")
    print(f"subset data rows:            {n_rows:,}")
    print(f"unique ids in subset:        {len(seen):,}")
    print(f"expected ids (from p/a):     {len(want):,}")
    print(f"COVERAGE  missing ids:       {len(missing):,}"
          + ("  <-- FAIL" if missing else "  OK"))
    print(f"SEQUENCES empty/NA seq rows: {n_empty_seq:,}"
          + ("  <-- FAIL" if n_empty_seq else "  OK"))
    print(f"INTEGRITY malformed rows:    {bad_rows:,}"
          + ("  <-- FAIL (truncated/partial write?)" if bad_rows else "  OK"))
    print(f"duplicate id rows:           {dup:,}  (informational)")
    if extra:
        print(f"ids in subset not in p/a:    {len(extra):,}  (informational)")
    print("===============================================")

    if missing:
        show = list(missing)[:10]
        print(f"\nfirst missing ids: {show}", file=sys.stderr)
    if args.write_missing and missing:
        with open(args.write_missing, "w") as fh:
            fh.write("\n".join(sorted(missing)) + "\n")
        print(f"missing ids written to {args.write_missing}")

    ok = not missing and not n_empty_seq and not bad_rows
    if ok:
        print("\nPASS: every p/a gene id is present, all sequences filled, no malformed rows.")
        sys.exit(0)
    else:
        print("\nFAIL: see report above. The run likely stopped early or the file is truncated.",
              file=sys.stderr)
        sys.exit(1)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pa", required=True, help="core p/a the subset was built from")
    p.add_argument("--subset", required=True, help="output gene_data.csv subset to validate")
    p.add_argument("--id-col", default="annotation_id", help="id column (default: annotation_id)")
    p.add_argument("--seq-col", default="dna_sequence", help="sequence column (default: dna_sequence)")
    p.add_argument("--write-missing", default=None,
                   help="optional path to write the list of missing ids (for a resume/re-run)")
    p.add_argument("--progress", action="store_true", help="print a line every 5M rows scanned")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
