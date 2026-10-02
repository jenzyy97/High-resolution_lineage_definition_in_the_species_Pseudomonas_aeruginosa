#!/usr/bin/env python3
"""
write_gene_fastas.py  --  pan-genome cleanup pipeline, step 3

Given a CLEANED core-genes gene_presence_absence.csv (step 2 output) and the
Panaroo merge `gene_data.csv`, write one nucleotide FASTA per gene family
(ready for per-gene alignment, e.g. with MAFFT). No genome downloads needed:
the DNA comes straight from gene_data.csv, keyed by the gene ids already sitting
in the p/a cells.

Join key: the ids in the p/a genome cells are Panaroo 'annotation_id' locus
tags (e.g. NAKOIJ_00001), matched against gene_data.csv's 'annotation_id'
column; the sequence is that row's 'dna_sequence'.

Each family FASTA holds one sequence per genome that carries the gene, with the
GENOME id (p/a column name / panaroo_query_id) as the FASTA header.

Assumes the first 3 p/a columns are metadata (Gene, Non-unique Gene name,
Annotation) and the rest are genome columns.

Example
-------
  python write_gene_fastas.py \
      --pa core_gene_presence_absence.filtered.csv \
      --gene-data gene_data.csv \
      --out-dir mafft_input_fasta
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

N_META_COLS = 3
META_COLS = ["Gene", "Non-unique Gene name", "Annotation"]
CELL_SPLIT = ";"  # paralog separator; cleaned p/a should have none, but be safe


def genome_cols(pa: pd.DataFrame) -> list:
    return [c for c in pa.columns if c not in META_COLS]


def cell_ids(cell) -> list:
    if pd.isna(cell) or str(cell).strip() == "":
        return []
    return [p.strip() for p in str(cell).split(CELL_SPLIT) if p.strip()]


def load_gene_data_seqs(path: str, id_col: str, seq_col: str, wanted_ids: set,
                        chunksize: int = 100_000, progress_every: int = 20) -> dict:
    """Stream gene_data.csv in chunks, keeping only the id->dna for ids we need.
    Streaming keeps memory bounded even though gene_data.csv is huge. Prints a
    heartbeat every `progress_every` chunks so a long run visibly progresses."""
    import time
    n_want = len(wanted_ids)
    print(f"Loading sequences from {path} (streaming) ...", flush=True)
    print(f"  (heartbeat every {progress_every} chunks of {chunksize:,} rows)", flush=True)
    seqs: dict = {}
    n_rows = 0
    n_chunks = 0
    t0 = time.time()
    reader = pd.read_csv(path, usecols=[id_col, seq_col], dtype=str, chunksize=chunksize)
    for chunk in reader:
        n_chunks += 1
        n_rows += len(chunk)
        hit = chunk[chunk[id_col].isin(wanted_ids)]
        for aid, dna in zip(hit[id_col], hit[seq_col]):
            if pd.notna(dna) and str(dna).strip():
                seqs[aid] = str(dna).strip()
        if n_chunks % progress_every == 0:
            elapsed = time.time() - t0
            rate = n_rows / elapsed if elapsed else 0
            print(f"  [{elapsed:6.0f}s] scanned {n_rows:>12,} rows | "
                  f"seqs found {len(seqs):>9,}/{n_want:,} | {rate:,.0f} rows/s", flush=True)
    print(f"  scanned {n_rows:,} gene_data rows in {time.time()-t0:.0f}s; found sequences for "
          f"{len(seqs):,}/{n_want:,} needed ids.")
    return seqs


def run(args: argparse.Namespace) -> None:
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Loading cleaned presence/absence file: {args.pa} ...", flush=True)
    pa = pd.read_csv(args.pa, dtype=str)
    gcols = genome_cols(pa)
    families = pa["Gene"].tolist()
    print(f"File loaded! {len(families)} families x {len(gcols)} genomes.")

    # collect every (gene_id, genome, family) referenced in the p/a cells.
    # Vectorised: melt the genome block to long form, drop blanks, split ';'
    # and explode -- avoids a Python loop over the ~tens of millions of cells.
    print("Collecting gene ids from p/a cells (vectorised) ...", flush=True)
    long = pa.melt(id_vars=["Gene"], value_vars=gcols,
                   var_name="genome", value_name="cell")
    long = long.dropna(subset=["cell"])
    long = long[long["cell"].astype(str).str.strip() != ""]
    long["gene_id"] = long["cell"].astype(str).str.split(CELL_SPLIT)
    long = long.explode("gene_id")
    long["gene_id"] = long["gene_id"].str.strip()
    long = long[long["gene_id"] != ""]
    id_pairs = list(zip(long["gene_id"], long["genome"], long["Gene"]))
    wanted = set(long["gene_id"].unique())
    print(f"  {len(wanted):,} unique gene ids referenced across {len(families):,} families.")

    seqs = load_gene_data_seqs(args.gene_data, args.id_col, args.seq_col, wanted,
                               chunksize=args.chunksize, progress_every=args.progress_every)

    missing_ids = wanted - set(seqs)
    if missing_ids:
        show = list(missing_ids)[:10]
        print(f"[warn] {len(missing_ids):,} referenced ids had no sequence in gene_data.csv "
              f"(first few: {show}). Those entries are skipped.", file=sys.stderr)

    # group sequences by family -> list of (genome, dna)
    per_family: dict = {fam: [] for fam in families}
    for gid, genome, fam in id_pairs:
        dna = seqs.get(gid)
        if dna:
            per_family[fam].append((genome, dna))

    print(f"Writing per-family FASTA files to {args.out_dir} ...", flush=True)
    n_written = 0
    n_done = 0
    empty = []
    total = len(families)
    for fam in families:
        n_done += 1
        entries = per_family[fam]
        if not entries:
            empty.append(fam)
            continue
        out_path = os.path.join(args.out_dir, f"{fam}.fasta")
        with open(out_path, "w") as fh:
            for genome, dna in entries:
                fh.write(f">{genome}\n{dna}\n")
        n_written += 1
        if n_written % args.fasta_progress_every == 0:
            print(f"  ... wrote {n_written:,} FASTAs ({n_done:,}/{total:,} families processed)",
                  flush=True)
        n_written += 1

    print(f"\nWrote {n_written} family FASTA(s) to {args.out_dir}")
    if empty:
        print(f"[warn] {len(empty)} families had no sequences and were skipped: "
              f"{empty[:10]}" + (" ..." if len(empty) > 10 else ""), file=sys.stderr)
    print("Done!")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pa", required=True, help="cleaned core gene_presence_absence.csv (step 2 output)")
    p.add_argument("--gene-data", required=True, help="Panaroo merge gene_data.csv")
    p.add_argument("--out-dir", required=True, help="directory to write per-family FASTA files")
    p.add_argument("--id-col", default="annotation_id",
                   help="column in gene_data.csv matching the p/a cell ids (default: annotation_id)")
    p.add_argument("--seq-col", default="dna_sequence",
                   help="column in gene_data.csv with the sequence to write (default: dna_sequence)")
    p.add_argument("--chunksize", type=int, default=100_000,
                   help="rows per streaming chunk (default: 100000)")
    p.add_argument("--progress-every", type=int, default=20,
                   help="print a heartbeat every N chunks while streaming gene_data (default: 20)")
    p.add_argument("--fasta-progress-every", type=int, default=500,
                   help="print a heartbeat every N FASTAs written (default: 500)")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
