#!/usr/bin/env python3
"""
filter_core_genome.py  --  pan-genome cleanup pipeline, step 2

Given a CORE-genes Panaroo gene_presence_absence.csv and thresholds read off the
step-1 distributions, produce a cleaned sub-presence/absence file:

  1. PARALOGS: drop gene families with a paralog (';') in >= --paralog-threshold %
     of genomes (genuine multi-copy families); then in the surviving families,
     blank (NA) any remaining cell that still contains a ';' (spurious paralogs,
     likely fragmentation / pipeline artifacts).
  2. PSEUDO: drop families with a 'pseudo' id in >= --pseudo-threshold % of
     genomes; then blank remaining cells containing 'pseudo'.
  3. REFOUND: left untouched (matches the original pipeline).
  4. Normalise empty strings to NA.
  5. Drop GENOMES (columns) that end up < --min-fill % filled after the blanking
     (too fragmented / low core-gene content).
  6. Write the cleaned sub-p/a.

The presence/absence file is loaded exactly once. Thresholds are percentages of
the number of genomes (matching the step-1 preview table).

Assumes the first 3 columns are metadata (Gene, Non-unique Gene name,
Annotation) and the rest are genome columns.

Example
-------
  python filter_core_genome.py \
      --pa core_gene_presence_absence.csv \
      --paralog-threshold 2 --pseudo-threshold 5 \
      --min-fill 80 \
      --out core_gene_presence_absence.filtered.csv
"""

from __future__ import annotations

import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

N_META_COLS = 3
META_COLS = ["Gene", "Non-unique Gene name", "Annotation"]


def load_pa(path: str) -> pd.DataFrame:
    print(f"Loading core presence/absence file: {path} ...", flush=True)
    pa = pd.read_csv(path, dtype=str)
    n_genomes = pa.shape[1] - N_META_COLS
    print(f"File loaded! {pa.shape[0]} gene families x {n_genomes} genomes.", flush=True)
    return pa


def genome_cols(pa: pd.DataFrame) -> list:
    return [c for c in pa.columns if c not in META_COLS]


# --------------------------------------------------------------------------
# per-family flag counts
# --------------------------------------------------------------------------
def paralog_counts(genomes: pd.DataFrame) -> pd.Series:
    """Per family: number of genomes whose cell contains a paralog (';')."""
    return genomes.apply(lambda col: col.astype(str).str.contains(";", na=False)).sum(axis=1)


def flag_counts(genomes: pd.DataFrame, marker: str) -> pd.Series:
    """Per family: number of genomes with >=1 id containing `marker`."""
    marker = marker.lower()

    def cell_flagged(cell) -> bool:
        if pd.isna(cell) or cell == "":
            return False
        return any(marker in part.lower() for part in str(cell).split(";"))

    return genomes.apply(lambda col: col.map(cell_flagged)).sum(axis=1)


# --------------------------------------------------------------------------
# filtering
# --------------------------------------------------------------------------
def drop_families_over_threshold(pa: pd.DataFrame, counts: pd.Series, min_g: int, label: str):
    """Drop family rows whose flag count >= min_g. Returns (kept_df, dropped_family_names)."""
    to_drop = counts >= min_g
    dropped_names = pa.loc[to_drop, "Gene"].tolist()
    kept = pa.loc[~to_drop]
    print(f"  {label}: dropped {len(dropped_names)} families (flag in >= {min_g} genomes); "
          f"{len(kept)} families remain.")
    return kept, dropped_names


def blank_cells(pa: pd.DataFrame, gcols: list, predicate, label: str):
    """Set to NA every genome cell for which predicate(cell)->True.
    Returns (new_df, n_cells_blanked)."""
    out = pa.copy()
    mask = out[gcols].apply(lambda col: col.map(predicate))
    out[gcols] = out[gcols].where(~mask, other=pd.NA)
    n = int(mask.values.sum())
    print(f"  {label}: blanked {n} cells.")
    return out, n


def cell_has_paralog(cell) -> bool:
    return (not pd.isna(cell)) and (";" in str(cell))


def make_marker_predicate(marker: str):
    marker = marker.lower()

    def pred(cell) -> bool:
        if pd.isna(cell) or cell == "":
            return False
        return any(marker in part.lower() for part in str(cell).split(";"))

    return pred


def plot_fill_before_after(pct_before: pd.Series, pct_after: pd.Series,
                           min_fill: float, out_path: str) -> None:
    """Two-panel bar of % filled families per genome (sorted), before and after
    the genome-fill drop. 'before' = after cell-blanking, all genomes; 'after' =
    genomes kept. The min-fill cutoff is drawn as a red line on the 'before' panel."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5), sharey=True)
    fig.suptitle("Filled cells per genome after cleaning", fontsize=13)

    for ax, pct, title in [
        (ax1, pct_before, f"Before genome drop (n={len(pct_before)})"),
        (ax2, pct_after, f"After genome drop (n={len(pct_after)})"),
    ]:
        vals = sorted(pct.values)
        ax.bar(range(len(vals)), vals, color="steelblue", edgecolor="none")
        mean = pct.mean() if len(pct) else 0
        ax.axhline(mean, color="black", linestyle="--", linewidth=1, label=f"Mean: {mean:.2f}%")
        ax.set_xlabel("Genome rank", fontsize=11)
        ax.set_title(title, fontsize=11)
        ax.legend(fontsize=9)
    # cutoff line only on the 'before' panel (that's where the drop is decided)
    ax1.axhline(min_fill, color="crimson", linestyle=":", linewidth=1.2,
                label=f"min-fill cutoff: {min_fill:.0f}%")
    ax1.legend(fontsize=9)
    ax1.set_ylabel("% gene families with valid entry", fontsize=11)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f".png saved to {out_path}")


def write_log(path: str, log: dict) -> None:
    """Write a human-readable, reproducible summary of the filtering."""
    lines = []
    lines.append("# filter_core_genome.py  --  filtering log")
    lines.append(f"generated: {log['timestamp']}")
    lines.append("")
    lines.append("## inputs")
    lines.append(f"input p/a:              {log['pa_path']}")
    lines.append(f"output p/a:             {log['out_path']}")
    lines.append(f"starting families:      {log['n_families_start']}")
    lines.append(f"starting genomes:       {log['n_genomes_start']}")
    lines.append("")
    lines.append("## thresholds")
    lines.append(f"paralog: >= {log['paralog_pct']}%  ({log['paralog_min']} genomes)")
    lines.append(f"pseudo:  >= {log['pseudo_pct']}%  ({log['pseudo_min']} genomes)")
    lines.append(f"refound: untouched")
    lines.append(f"min genome fill: >= {log['min_fill_pct']}%")
    lines.append("")
    lines.append("## paralog step")
    lines.append(f"families dropped (genuine paralogs): {len(log['paralog_families_dropped'])}")
    lines.append(f"spurious paralog cells blanked:      {log['paralog_cells_blanked']}")
    lines.append("")
    lines.append("## pseudo step")
    lines.append(f"families dropped (genuine pseudo):   {len(log['pseudo_families_dropped'])}")
    lines.append(f"pseudo cells blanked:                {log['pseudo_cells_blanked']}")
    lines.append("")
    lines.append("## genome fill step")
    lines.append(f"genomes kept:    {log['n_genomes_kept']}")
    lines.append(f"genomes dropped: {len(log['genomes_dropped'])}")
    lines.append("")
    lines.append("## result")
    lines.append(f"final families: {log['n_families_final']}")
    lines.append(f"final genomes:  {log['n_genomes_final']}")
    lines.append("")
    lines.append("## dropped family names (paralog step)")
    lines.extend(log["paralog_families_dropped"] or ["(none)"])
    lines.append("")
    lines.append("## dropped family names (pseudo step)")
    lines.extend(log["pseudo_families_dropped"] or ["(none)"])
    lines.append("")
    lines.append("## dropped genome names (fill < threshold)")
    lines.extend(log["genomes_dropped"] or ["(none)"])
    lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    print(f"Log written -> {path}")


def run(args: argparse.Namespace) -> None:
    import datetime

    pa = load_pa(args.pa)
    gcols = genome_cols(pa)
    n_genomes = len(gcols)

    par_min = int(n_genomes * args.paralog_threshold / 100)
    pse_min = int(n_genomes * args.pseudo_threshold / 100)
    print(f"\nThresholds: paralog >= {args.paralog_threshold}% ({par_min} genomes), "
          f"pseudo >= {args.pseudo_threshold}% ({pse_min} genomes). "
          f"refound: untouched.\n")

    log = {
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
        "pa_path": args.pa, "out_path": args.out,
        "n_families_start": pa.shape[0], "n_genomes_start": n_genomes,
        "paralog_pct": args.paralog_threshold, "paralog_min": par_min,
        "pseudo_pct": args.pseudo_threshold, "pseudo_min": pse_min,
        "min_fill_pct": args.min_fill,
    }

    # 1. paralogs: drop genuine families, then blank spurious ';' cells
    print("Paralog step:")
    par_c = paralog_counts(pa[gcols])
    pa, log["paralog_families_dropped"] = drop_families_over_threshold(
        pa, par_c, par_min, "paralog family drop")
    gcols = genome_cols(pa)
    pa, log["paralog_cells_blanked"] = blank_cells(
        pa, gcols, cell_has_paralog, "spurious paralog cells")

    # 2. pseudo: drop genuine families, then blank remaining pseudo cells
    print("Pseudo step:")
    pse_c = flag_counts(pa[gcols], "pseudo")
    pa, log["pseudo_families_dropped"] = drop_families_over_threshold(
        pa, pse_c, pse_min, "pseudo family drop")
    pa, log["pseudo_cells_blanked"] = blank_cells(
        pa, gcols, make_marker_predicate("pseudo"), "remaining pseudo cells")

    # 3. refound: untouched (matches original)

    # 4. normalise empty strings to NA
    pa[gcols] = pa[gcols].replace("", pd.NA)

    # 5. drop genomes (columns) below the fill threshold
    print("\nGenome fill step:")
    frac_filled = pa[gcols].notna().mean()          # fraction 0..1, per genome (after blanking)
    pct_before = frac_filled * 100
    keep = frac_filled[frac_filled >= args.min_fill / 100].index.tolist()
    drop = frac_filled[frac_filled < args.min_fill / 100].index.tolist()
    pct_after = pct_before[keep]
    print(f"  min-fill = {args.min_fill}% -> keeping {len(keep)} genomes, dropping {len(drop)}.")
    print(f"  %% filled per genome BEFORE drop (min/mean/max): "
          f"{pct_before.min():.2f}% / {pct_before.mean():.2f}% / {pct_before.max():.2f}%")
    if len(pct_after):
        print(f"  %% filled per genome AFTER  drop (min/mean/max): "
              f"{pct_after.min():.2f}% / {pct_after.mean():.2f}% / {pct_after.max():.2f}%")
    if drop:
        show = drop[:10]
        print(f"  dropped genomes (first {len(show)}): {show}" + (" ..." if len(drop) > 10 else ""))
    pa = pa[META_COLS + keep]
    log["genomes_dropped"] = drop
    log["n_genomes_kept"] = len(keep)

    # 6. write outputs
    pa.to_csv(args.out, index=False)
    log["n_families_final"] = pa.shape[0]
    log["n_genomes_final"] = len(keep)
    print(f"\nCleaned sub-p/a written: {pa.shape[0]} families x {len(keep)} genomes -> {args.out}")

    stem = args.out.rsplit(".", 1)[0]
    log_path = args.log if args.log else stem + ".filter_log.txt"
    write_log(log_path, log)

    fig_path = args.fill_fig if args.fill_fig else stem + ".fill_per_genome.png"
    plot_fill_before_after(pct_before, pct_after, args.min_fill, fig_path)
    print("Done!")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pa", required=True, help="core-genes Panaroo gene_presence_absence.csv")
    p.add_argument("--out", required=True, help="output path for the cleaned sub-p/a CSV")
    p.add_argument("--paralog-threshold", type=float, required=True,
                   help="drop families with a paralog in >= this %% of genomes (genuine paralogs); "
                        "remaining ';' cells are blanked")
    p.add_argument("--pseudo-threshold", type=float, required=True,
                   help="drop families with a pseudo id in >= this %% of genomes; remaining pseudo "
                        "cells are blanked")
    p.add_argument("--min-fill", type=float, default=80.0,
                   help="drop genomes (columns) filled in < this %% of families after blanking "
                        "(default: 80)")
    p.add_argument("--log", default=None,
                   help="path for the filtering log (default: <out>.filter_log.txt next to --out)")
    p.add_argument("--fill-fig", default=None,
                   help="path for the before/after fill-per-genome QC plot "
                        "(default: <out>.fill_per_genome.png next to --out)")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
