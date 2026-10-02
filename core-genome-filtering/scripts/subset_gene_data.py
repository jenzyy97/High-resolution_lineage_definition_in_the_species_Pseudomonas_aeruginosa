#!/usr/bin/env python3
"""
subset_gene_data.py  --  pan-genome cleanup pipeline (companion)

Write a subset of the Panaroo merge `gene_data.csv` containing only the rows for
the gene families in a given (core) p/a -- an alternative, much smaller
gene_data.csv to ship alongside the core pan-genome.

Every gene id appearing in the p/a's genome cells is gathered, and every
gene_data.csv row whose id column is in that set is kept.

  --pa        the core p/a (e.g. the 98%-presence core gene_presence_absence.csv)
  --gene-data merge gene_data.csv (streamed) -> rows kept where id in the set

Assumes the first 3 p/a columns are metadata (Gene, Non-unique Gene name,
Annotation) and the rest are genome columns; gene_data ids are matched on
'annotation_id' by default.

Example
-------
  python subset_gene_data.py \
      --pa core_gene_presence_absence.csv \
      --gene-data gene_data.csv \
      --out gene_data.core.csv
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

META_COLS = ["Gene", "Non-unique Gene name", "Annotation"]
CELL_SPLIT = ";"


def genome_cols(pa: pd.DataFrame) -> list:
    return [c for c in pa.columns if c not in META_COLS]


def cell_ids(cell) -> list:
    if pd.isna(cell) or str(cell).strip() == "":
        return []
    return [p.strip() for p in str(cell).split(CELL_SPLIT) if p.strip()]


def all_cell_ids(pa: pd.DataFrame) -> set:
    """Every gene id appearing in the p/a's genome cells.

    Vectorised: flatten all genome columns to a 1-D Series, drop blanks, then
    split on ';' in bulk and explode. Avoids a Python loop over the ~tens of
    millions of individual cells (which is minutes vs. seconds at this scale)."""
    gcols = genome_cols(pa)
    # stack every genome column into one long Series of cell strings
    flat = pa[gcols].to_numpy().ravel()
    s = pd.Series(flat, dtype="object").dropna()
    s = s[s.astype(str).str.strip() != ""]
    if s.empty:
        return set()
    # split paralog cells and flatten; strip whitespace; drop empties
    parts = s.astype(str).str.split(CELL_SPLIT).explode()
    parts = parts.str.strip()
    parts = parts[parts != ""]
    return set(parts.unique())


def run(args: argparse.Namespace) -> None:
    print(f"Loading core p/a: {args.pa} ...", flush=True)
    pa = pd.read_csv(args.pa, dtype=str)
    print(f"  {pa.shape[0]} families x {len(genome_cols(pa))} genomes.")

    print("Gathering all gene ids from the p/a cells ...", flush=True)
    wanted = all_cell_ids(pa)
    print(f"  {len(wanted)} unique gene ids across {pa.shape[0]} families.")

    n_want = len(wanted)
    print(f"Streaming {args.gene_data} and keeping core rows ...", flush=True)
    print(f"  (heartbeat every {args.progress_every} chunks of {args.chunksize:,} rows)",
          flush=True)
    import time
    t0 = time.time()
    n_in, n_out, wrote_header = 0, 0, False
    found_ids: set = set()
    n_chunks = 0
    reader = pd.read_csv(args.gene_data, dtype=str, chunksize=args.chunksize)
    for chunk in reader:
        n_chunks += 1
        n_in += len(chunk)
        if args.id_col not in chunk.columns:
            sys.exit(f"id column '{args.id_col}' not in gene_data.csv "
                     f"(columns: {list(chunk.columns)})")
        keep = chunk[chunk[args.id_col].isin(wanted)]
        if len(keep):
            keep.to_csv(args.out, mode="a" if wrote_header else "w",
                        header=not wrote_header, index=False)
            wrote_header = True
            n_out += len(keep)
            found_ids.update(keep[args.id_col].dropna())
        if n_chunks % args.progress_every == 0:
            elapsed = time.time() - t0
            rate = n_in / elapsed if elapsed else 0
            print(f"  [{elapsed:6.0f}s] scanned {n_in:>12,} rows | kept {n_out:>11,} | "
                  f"ids found {len(found_ids):>9,}/{n_want:,} | {rate:,.0f} rows/s",
                  flush=True)

    if not wrote_header:
        # nothing matched: still write an empty file with the header
        pd.read_csv(args.gene_data, nrows=0).to_csv(args.out, index=False)

    elapsed = time.time() - t0
    print(f"\nScanned {n_in:,} gene_data rows in {elapsed:.0f}s; "
          f"wrote {n_out:,} core rows -> {args.out}")

    # coverage from the ids we actually matched during streaming (no re-read needed)
    missing = wanted - found_ids
    print(f"Coverage: {len(found_ids):,}/{n_want:,} core ids matched a gene_data row.")
    if missing:
        print(f"[warn] {len(missing):,} core ids had NO gene_data row "
              f"(first few: {list(missing)[:5]}).", file=sys.stderr)
        if args.write_missing:
            with open(args.write_missing, "w") as fh:
                fh.write("\n".join(sorted(missing)) + "\n")
            print(f"       missing ids written to {args.write_missing}")
    print("Done!")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pa", required=True, help="core p/a (source of gene ids to keep)")
    p.add_argument("--gene-data", required=True, help="Panaroo merge gene_data.csv (streamed)")
    p.add_argument("--out", required=True, help="output path for the subset gene_data.csv")
    p.add_argument("--id-col", default="annotation_id",
                   help="gene_data.csv column matching the p/a cell ids (default: annotation_id)")
    p.add_argument("--chunksize", type=int, default=100_000,
                   help="rows per streaming chunk (default: 100000)")
    p.add_argument("--progress-every", type=int, default=20,
                   help="print a heartbeat every N chunks (default: 20)")
    p.add_argument("--write-missing", default=None,
                   help="if any core ids have no gene_data row, write them to this path")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
