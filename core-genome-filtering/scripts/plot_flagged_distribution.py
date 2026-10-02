#!/usr/bin/env python3
"""
plot_flagged_distribution.py  --  pan-genome cleanup pipeline, step 1

Given a CORE-genes Panaroo gene_presence_absence.csv (already subset to the
gene families deemed core), plot the distribution of:
  (1) paralog-containing genomes per gene family  (cells holding >1 id, ';')
  (2) pseudo- and refound-flagged genomes per gene family

These distributions help you choose thresholds for the next step (filtering).
The presence/absence file is loaded exactly once.

Optionally restrict to a SUBSET of genomes (--subset) before any counting, so
the distributions reflect only the genomes you care about. Ids in the subset
file may be exact panaroo_query_ids, GCF/GCA accessions, or strain names; they
are resolved directly against the p/a genome column names (the accession core
and strain name are embedded there, so no metadata table is needed).

Assumes the first 3 columns are metadata (Gene, Non-unique Gene name,
Annotation) and the rest are genome columns.

Example
-------
  # whole core pan-genome
  python plot_flagged_distribution.py --pa core_gene_presence_absence.csv --out-dir ./qc

  # restricted to a subset of genomes
  python plot_flagged_distribution.py --pa core_gene_presence_absence.csv --out-dir ./qc \
      --subset my_genomes.txt
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

N_META_COLS = 3  # leading metadata columns before genome columns
GCF_GCA_RE = re.compile(r"GC[AF]_?(\d+)", re.IGNORECASE)


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def load_pa(path: str) -> pd.DataFrame:
    print(f"Loading core presence/absence file: {path} ...", flush=True)
    pa = pd.read_csv(path, dtype=str)
    n_genomes = pa.shape[1] - N_META_COLS
    print(f"File loaded! {pa.shape[0]} gene families x {n_genomes} genomes.", flush=True)
    return pa


def genome_columns(pa: pd.DataFrame) -> list:
    return list(pa.columns[N_META_COLS:])


# --------------------------------------------------------------------------
# genome subset resolution (exact id / GCF-GCA accession core / strain name)
# --------------------------------------------------------------------------
def _normalise_strain(s: str) -> str:
    return re.sub(r"[\s_.\-]+", "", str(s).lower())


def resolve_one(query: str, gcols: list) -> str:
    """Resolve a single id to an exact p/a genome column (or raise)."""
    if query in gcols:
        return query
    m = GCF_GCA_RE.search(query)
    if m:
        core = m.group(1)
        hits = [c for c in gcols if core in c]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            raise ValueError(f"accession core '{core}' matched multiple genomes: {hits[:10]}")
    q = _normalise_strain(query)
    if q:
        hits = [c for c in gcols if q in _normalise_strain(c)]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            raise ValueError(f"strain query '{query}' matched multiple genomes: {hits[:10]}")
    raise ValueError(f"could not resolve '{query}' to any genome column")


def read_subset_file(path: str) -> list:
    with open(path) as fh:
        ids = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    if not ids:
        sys.exit(f"subset file '{path}' contained no ids.")
    return ids


def resolve_subset(pa: pd.DataFrame, subset_ids: list) -> list:
    """Map each subset id to a p/a genome column; report failures, keep order/unique."""
    gcols = genome_columns(pa)
    resolved, unresolved, seen = [], [], set()
    for q in subset_ids:
        try:
            col = resolve_one(q, gcols)
        except ValueError as e:
            unresolved.append((q, str(e)))
            continue
        if col not in seen:
            seen.add(col)
            resolved.append(col)
    print(f"Subset: {len(resolved)}/{len(subset_ids)} ids resolved to genome columns.")
    if unresolved:
        print(f"[warn] {len(unresolved)} id(s) could not be resolved:", file=sys.stderr)
        for q, msg in unresolved:
            print(f"        - {q}: {msg}", file=sys.stderr)
    if not resolved:
        sys.exit("No subset ids resolved to genome columns; nothing to plot.")
    return resolved


# --------------------------------------------------------------------------
# flag counting
# --------------------------------------------------------------------------
def paralog_genome_counts(genomes: pd.DataFrame) -> pd.Series:
    """Per gene family (row): number of genomes whose cell holds a paralog (';')."""
    has_paralog = genomes.apply(lambda col: col.astype(str).str.contains(";", na=False))
    return has_paralog.sum(axis=1)


def flag_genome_counts(genomes: pd.DataFrame, marker: str) -> pd.Series:
    """Per gene family (row): number of genomes with >=1 id containing `marker`."""
    marker = marker.lower()

    def cell_flagged(cell) -> bool:
        if pd.isna(cell) or cell == "":
            return False
        return any(marker in part.lower() for part in str(cell).split(";"))

    flagged = genomes.apply(lambda col: col.map(cell_flagged))
    return flagged.sum(axis=1)


# --------------------------------------------------------------------------
# auto split: boundary = first genome-count bin (>=1) with zero families
# --------------------------------------------------------------------------
def auto_split(counts: pd.Series) -> int | None:
    """Given per-family paralog counts, return the left edge of the valley that
    ends the first (spurious) peak: the peak's tallest genome-count bin, then the
    first empty bin *after* that peak. Walking only after the peak has been
    reached makes this robust to small gaps inside the rising edge. Returns None
    if the non-zero counts run contiguously to the maximum with no gap past the
    peak."""
    dist = counts[counts > 0].value_counts().sort_index()  # index=genome-count, val=#families
    if dist.empty:
        return None
    present = set(int(k) for k in dist.index)
    peak_count = int(dist.idxmax())  # genome-count value of the tallest bar
    hi = max(present)
    for k in range(peak_count + 1, hi + 1):
        if k not in present:
            return k
    return None  # no empty bin after the peak


def _add_pct_axis(ax, n_genomes: int):
    axx = ax.twiny()
    axx.set_xlim(ax.get_xlim())
    x_min, x_max = ax.get_xlim()
    raw_ticks = [p / 100 * n_genomes for p in range(0, 101, 25)]
    visible = [(r, p) for r, p in zip(raw_ticks, range(0, 101, 25)) if x_min <= r <= x_max]
    if visible:
        raw, pct = zip(*visible)
        axx.set_xticks(raw)
        axx.set_xticklabels([f"{p}%" for p in pct])
    axx.set_xlabel("Percentage of genomes", fontsize=11)
    return axx


def plot_paralog_distribution(counts: pd.Series, n_genomes: int, split: int, out_path: str):
    """Two-panel low/high split at `split`; left = first (spurious) peak,
    right = valley + genuine-paralog tail. A vertical line marks the split."""
    counts = counts[counts > 0]
    dist = counts.value_counts().sort_index()
    low = dist[dist.index < split]
    high = dist[dist.index >= split]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))
    fig.suptitle("Distribution of paralog-containing genomes per gene family", fontsize=13)

    ax1.bar(low.index, low.values, width=2.0, color="steelblue", edgecolor="steelblue", linewidth=0.5)
    ax1.set_xlabel("Number of genomes with paralogs", fontsize=11)
    ax1.set_ylabel("Number of gene families", fontsize=11)
    ax1.set_title(f"First peak (< {split} genomes, n={int(low.values.sum())})", fontsize=11)
    ax1.axvline(split, color="crimson", linestyle="--", linewidth=1, label=f"split = {split}")
    ax1.legend(fontsize=9)
    _add_pct_axis(ax1, n_genomes)

    ax2.bar(high.index, high.values, width=2.0, color="darkorange", edgecolor="darkorange", linewidth=0.5)
    ax2.set_xlabel("Number of genomes with paralogs", fontsize=11)
    ax2.set_ylabel("Number of gene families", fontsize=11)
    pct = split / n_genomes * 100 if n_genomes else 0
    ax2.set_title(f"Genuine tail (>= {split}, ~{pct:.1f}% genomes, n={int(high.values.sum())})", fontsize=11)
    _add_pct_axis(ax2, n_genomes)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f".png saved to {out_path}")


def plot_flag_distribution(pseudo: pd.Series, refound: pd.Series, n_genomes: int,
                           out_path: str, bin_width: int = 1):
    """One histogram per flag, one bar per genome-count value (bin width 1).
    Counts can span 0..n_genomes but the mass sits at the low end, so each panel
    auto-zooms its x-axis to the 99th percentile of its own flagged counts -- this
    keeps the thin bars visible instead of squashing them against a long tail.
    A note in the title says how many families fall beyond the zoom."""
    bin_width = max(1, int(bin_width))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))
    fig.suptitle("Distribution of flagged genomes per gene family", fontsize=13)

    for ax, counts, label, color in [
        (ax1, pseudo, "pseudo", "steelblue"),
        (ax2, refound, "refound", "darkorange"),
    ]:
        vals = counts[counts > 0].to_numpy()
        n_fam = len(vals)
        beyond = 0
        if n_fam:
            hi = int(vals.max())
            # auto-zoom: 99th percentile, but never below a small floor so a very
            # tight distribution still shows a sensible range
            p99 = int(np.ceil(np.percentile(vals, 99)))
            xmax = max(10, min(hi, p99))
            beyond = int((vals > xmax).sum())
            edges = list(range(1, xmax + bin_width + 1, bin_width))
            if len(edges) < 2:
                edges = [1, 1 + bin_width]
            ax.hist(vals, bins=edges, color=color, edgecolor=color, linewidth=0.3)
            ax.set_xlim(0, xmax)
        ax.set_xlabel("Number of genomes with flag", fontsize=11)
        ax.set_ylabel("Number of gene families", fontsize=11)
        title = f"'{label}' flag (n={n_fam} families with >=1 flag)"
        if beyond:
            title += f"\n(+{beyond} families beyond x-axis zoom)"
        ax.set_title(title, fontsize=11)
        _add_pct_axis(ax, n_genomes)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f".png saved to {out_path}")


def print_threshold_helper(name: str, counts: pd.Series, n_genomes: int):
    """Preview of how many families each candidate threshold would drop. The
    threshold itself is chosen and applied in the next (filtering) step."""
    counts = counts[counts > 0]
    print(f"\n{name}: {len(counts)} families have >=1 flagged genome.")
    print(f"  If you set the threshold to X%, families flagged in >= X% of "
          f"{n_genomes} genomes would be dropped:")
    for pct in (0.1, 0.5, 1, 2, 5, 10, 25, 50):
        min_g = int(n_genomes * pct / 100)
        n_drop = int((counts >= min_g).sum())
        print(f"    >= {pct:>4}%  (>= {min_g:>5} genomes):  {n_drop} families dropped")


def run(args: argparse.Namespace) -> None:
    os.makedirs(args.out_dir, exist_ok=True)
    pa = load_pa(args.pa)

    # optional genome subset (applied before any counting)
    if args.subset:
        subset_ids = read_subset_file(args.subset)
        cols = resolve_subset(pa, subset_ids)
        genomes = pa[cols]
    else:
        genomes = pa[genome_columns(pa)]
    n_genomes = genomes.shape[1]
    print(f"Using {n_genomes} genomes for the distributions.\n")

    print("Counting paralog cells per family ...", flush=True)
    par_counts = paralog_genome_counts(genomes)
    print("Counting pseudo/refound cells per family ...", flush=True)
    pseudo_counts = flag_genome_counts(genomes, "pseudo")
    refound_counts = flag_genome_counts(genomes, "refound")

    # decide split
    if args.paralog_split is not None:
        split = args.paralog_split
        print(f"\nUsing manual paralog split = {split} genomes.")
    else:
        split = auto_split(par_counts)
        if split is None:
            split = max(1, int(n_genomes * 0.1))
            print(f"\n[warn] no empty bin found in the paralog distribution; "
                  f"falling back to split = {split} (10% of genomes). "
                  f"Override with --paralog-split.", file=sys.stderr)
        else:
            pct = split / n_genomes * 100
            print(f"\nAuto paralog split = {split} genomes ({pct:.2f}% of genomes) "
                  f"-- first empty bin after the low peak.")
            print(f"Suggested paralog threshold for filtering: {split} genomes "
                  f"({pct:.2f}%) -- families with paralogs in >= {split} genomes look genuine; "
                  f"below that looks spurious.")

    plot_paralog_distribution(
        par_counts, n_genomes, split,
        os.path.join(args.out_dir, "paralog_distribution.png"),
    )
    plot_flag_distribution(
        pseudo_counts, refound_counts, n_genomes,
        os.path.join(args.out_dir, "pseudo_refound_distribution.png"),
        bin_width=args.flag_bin_width,
    )

    print_threshold_helper("paralogs", par_counts, n_genomes)
    print_threshold_helper("pseudo", pseudo_counts, n_genomes)
    print_threshold_helper("refound", refound_counts, n_genomes)

    print("\nDone! Use these distributions to pick thresholds for the filtering step.")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pa", required=True, help="core-genes Panaroo gene_presence_absence.csv")
    p.add_argument("--out-dir", default=".", help="directory for the output PNGs (default: current dir)")
    p.add_argument("--paralog-split", type=int, default=None,
                   help="manual low/high panel split (genomes). Default: auto-detect the first "
                        "empty bin after the low peak.")
    p.add_argument("--flag-bin-width", type=int, default=1,
                   help="bin width (in genomes) for the pseudo/refound histograms (default: 1, "
                        "one bar per genome-count value). The x-axis auto-zooms to the 99th "
                        "percentile of each flag's counts so the low-end bars stay visible.")
    g = p.add_argument_group("optional genome subset")
    g.add_argument("--subset", help="file of genome ids (one per line: exact panaroo_query_id, "
                                    "GCF/GCA accession, or strain name). Ids are resolved against "
                                    "the p/a column names. Restricts the analysis to these "
                                    "genomes. Default: all genomes.")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
