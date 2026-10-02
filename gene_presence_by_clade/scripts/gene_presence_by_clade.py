#!/usr/bin/env python3
"""
gene_presence_by_clade.py

Given a Panaroo gene_presence_absence.csv and a clade assignment table,
plot how many genomes per clade carry a chosen gene, and optionally write
a CSV of the genomes that have it (with clade + specific gene id).

The gene can be specified either as:
  * a Panaroo gene family  (a value in the 'Gene' column), with --gene-family, or
  * a specific gene id found in one genome's column, with a genome + --gene-id
    (the script resolves which family that id belongs to). A leading prefix such
    as 'cds-' is matched tolerantly, so 'cds-WP_012076956.1' and the bare
    'WP_012076956.1' both work.

The genome for the --gene-id path can be given as:
  * --genome  : an EXACT panaroo genome id (matches a column in the p/a file,
                i.e. the 'panaroo_query_id' value), or
  * --query   : a looser identifier the script resolves to a panaroo genome id:
                a GCF/GCA accession (the numeric core, e.g. 000017205 from
                GCF_000017205.1, is matched against the p/a column names), or a
                strain name (normalised match against the column names).

Not sure what a genome's gene ids look like? Inspect first with --list-ids:
  python gene_presence_by_clade.py --pa gene_presence_absence.csv \
      --query GCF_000017205.1 --list-ids

Example
-------
  # by panaroo gene family
  python gene_presence_by_clade.py \
      --pa gene_presence_absence.csv \
      --clades All_Pae_metadata.csv \
      --panaroo-id-col panaroo_query_id --clade-col clade \
      --gene-family group_36431 \
      --out-fig group_36431_by_clade.png --out-csv group_36431_genomes.csv

  # by a specific gene id, exact panaroo genome id
  python gene_presence_by_clade.py \
      --pa gene_presence_absence.csv \
      --clades Pae_popPUNKcluster_representatives_metadata.csv \
      --panaroo-id-col panaroo_query_id --clade-col clade \
      --genome Pseudomonas_aeruginosa_PA7_119_GCF_000017205.1 \
      --gene-id cds-WP_012076956.1 --out-fig out.png

  # by a specific gene id, genome found via GCF/GCA accession or strain
  python gene_presence_by_clade.py \
      --pa gene_presence_absence.csv \
      --clades All_Pae_metadata.csv \
      --panaroo-id-col panaroo_query_id --clade-col clade \
      --query GCF_000017205.1 --gene-id cds-WP_012076956.1 --out-fig out.png
"""

from __future__ import annotations

import argparse
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Columns in a standard Panaroo/roary gene_presence_absence.csv that are
# metadata, i.e. NOT genome columns.
DEFAULT_METADATA_COLS = [
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
    "No. sequences", "Avg sequences per isolate", "Genome Fragment",
    "Order within Fragment", "Accessory Fragment",
    "Accessory Order with Fragment", "QC", "Min group size nuc",
    "Max group size nuc", "Avg group size nuc",
]

CELL_SPLIT = re.compile(r"[;\t]")

# Clade labels dropped from the FIGURE by default (they stay in the output CSV
# and the printed summary). These are non-biological "clades" from the metadata.
DEFAULT_PLOT_EXCLUDE = [
    "Genome removed from analysis due to low core gene content (high fragmentation)",
]

# Panaroo/prokka often prefix ids (e.g. 'cds-WP_012074645.1' where the real
# protein id is 'WP_012074645.1'). Strip a leading '<word>-' so a user pasting
# the bare id from a GFF still matches.
ID_PREFIX_RE = re.compile(r"^[A-Za-z]+-")


def normalise_gene_id(gid: str) -> str:
    """Strip a leading '<prefix>-' (e.g. 'cds-') for tolerant matching."""
    return ID_PREFIX_RE.sub("", str(gid).strip())


def genome_columns(pa: pd.DataFrame, metadata_cols=None) -> list:
    metadata_cols = metadata_cols or DEFAULT_METADATA_COLS
    return [c for c in pa.columns if c not in metadata_cols]


# Substrings that mark a Panaroo-synthesised id rather than a real annotated
# gene id (re-found genes, pseudogenes, length-placeholder stubs). Matched
# case-insensitively against each id.
NON_REAL_ID_MARKERS = ("refound", "pseudo", "len_")


def cell_ids(cell) -> list:
    """Split a presence/absence cell into individual gene ids (handles NaN,
    paralogs separated by ';' or tab)."""
    if pd.isna(cell):
        return []
    return [p.strip() for p in CELL_SPLIT.split(str(cell)) if p.strip()]


def is_real_gene_id(gid: str) -> bool:
    """True if gid looks like a real annotated gene id, i.e. not a Panaroo
    'refound'/'pseudo'/length-placeholder id and not blank."""
    if not gid or not str(gid).strip():
        return False
    low = str(gid).lower()
    return not any(marker in low for marker in NON_REAL_ID_MARKERS)


def id_examples(pa: pd.DataFrame, genome: str, n: int = 20, real_only: bool = True) -> list:
    """Up to n example gene ids from a genome's column, to show id format.
    By default skips NaN and Panaroo synthetic ids (refound/pseudo/len_)."""
    if genome not in pa.columns:
        raise KeyError(f"'{genome}' is not a column in the presence/absence file.")
    out = []
    for cell in pa[genome].dropna():
        for gid in cell_ids(cell):
            if real_only and not is_real_gene_id(gid):
                continue
            out.append(gid)
            if len(out) >= n:
                return out
    return out


def print_id_examples(pa: pd.DataFrame, genome: str, n: int = 20) -> None:
    """Inspection helper: print up to n real gene ids from a genome's column
    so the user can see the id format before choosing a --gene-id."""
    if genome not in pa.columns:
        gcols = genome_columns(pa)
        raise KeyError(
            f"'{genome}' is not a genome column in the presence/absence file.\n"
            f"  Example genome columns: {gcols[:3]}"
        )
    real = id_examples(pa, genome, n=n, real_only=True)
    if not real:
        print(f"No real (non-refound/pseudo) gene ids found in genome '{genome}'.")
        return
    print(f"Up to {n} example gene ids in genome '{genome}' "
          f"(refound/pseudo/len_ and blanks excluded):")
    for gid in real:
        stripped = normalise_gene_id(gid)
        note = f"   (bare id: {stripped})" if stripped != gid else ""
        print(f"  {gid}{note}")


def _id_signature(gid: str) -> str:
    m = re.match(r"^([A-Za-z_.\-]*)", gid)
    prefix = m.group(1) if m else ""
    return f"prefix={prefix!r}|parts={len(gid.split('_'))}"


def check_id_format(pa: pd.DataFrame, genome: str, gene_id: str, n: int = 5) -> bool:
    """Warn (don't fail) if gene_id looks structurally unlike the genome's
    other ids. Returns True if it looks consistent."""
    examples = id_examples(pa, genome, n=n)
    if not examples:
        print(f"[warn] no example ids in genome '{genome}' to compare against.", file=sys.stderr)
        return False
    ok = _id_signature(gene_id) in {_id_signature(e) for e in examples}
    if not ok:
        print(
            f"[warn] gene id '{gene_id}' doesn't match the format of ids in "
            f"genome '{genome}'.\n       examples: {examples}",
            file=sys.stderr,
        )
    return ok


def resolve_family(pa: pd.DataFrame, genome: str, gene_id: str) -> str:
    """Find which family ('Gene') a specific gene id belongs to in a genome.

    Matching tries an exact match first, then a prefix-tolerant match so a
    bare protein id (e.g. 'WP_012074645.1') resolves against a prefixed cell
    value (e.g. 'cds-WP_012074645.1') and vice versa.
    """
    if "Gene" not in pa.columns:
        raise KeyError("presence/absence file has no 'Gene' column.")
    if genome not in pa.columns:
        raise KeyError(f"'{genome}' is not a column in the presence/absence file.")

    # 1. exact match
    mask = pa[genome].apply(lambda cell: gene_id in cell_ids(cell))
    matches = pa.loc[mask, "Gene"].unique().tolist()

    # 2. prefix-tolerant fallback (strip 'cds-' etc. from both sides)
    if not matches:
        q_norm = normalise_gene_id(gene_id)
        mask = pa[genome].apply(
            lambda cell: q_norm in [normalise_gene_id(i) for i in cell_ids(cell)]
        )
        matches = pa.loc[mask, "Gene"].unique().tolist()
        if matches:
            print(
                f"[info] '{gene_id}' matched via prefix-tolerant lookup "
                f"(compared as '{q_norm}').",
                file=sys.stderr,
            )

    if not matches:
        raise ValueError(
            f"gene id '{gene_id}' not found in genome '{genome}'. "
            f"Example ids there: {id_examples(pa, genome)}"
        )
    if len(matches) > 1:
        print(f"[warn] '{gene_id}' matched multiple families {matches}; using first.", file=sys.stderr)
    return matches[0]


GCF_GCA_RE = re.compile(r"GC[AF]_?(\d+)", re.IGNORECASE)


def _normalise_strain(s: str) -> str:
    """Lowercase and strip separators so 'Pseudomonas aeruginosa PA7' and
    'pseudomonas_aeruginosa-pa7' compare equal."""
    return re.sub(r"[\s_.\-]+", "", str(s).lower())


def resolve_genome(pa: pd.DataFrame, query: str, metadata_cols=None) -> str:
    """Resolve a loose identifier to an exact panaroo genome id (a column in
    the p/a file, i.e. a 'panaroo_query_id' value).

    Match order:
      1. exact column match
      2. GCF/GCA accession -- the numeric core (e.g. 000017205 from
         GCF_000017205.1) is matched against the genome column names
      3. strain name -- normalised (case/space/_/-/. removed) substring match
         against the genome column names

    Raises if nothing matches or the match is ambiguous.
    """
    gcols = genome_columns(pa, metadata_cols)

    # 1. exact
    if query in gcols:
        return query

    # 2. GCF/GCA accession -> numeric core
    m = GCF_GCA_RE.search(query)
    if m:
        core = m.group(1)  # digits only, e.g. '000017205'
        hits = [c for c in gcols if core in c]
        if len(hits) == 1:
            print(f"Resolved accession '{query}' (core {core}) -> genome '{hits[0]}'")
            return hits[0]
        if len(hits) > 1:
            raise ValueError(
                f"accession core '{core}' matched multiple genomes: {hits[:10]}"
                + (" ..." if len(hits) > 10 else "")
            )
        # fall through to strain matching if the accession core matched nothing

    # 3. strain name -> normalised substring match
    q_norm = _normalise_strain(query)
    if q_norm:
        hits = [c for c in gcols if q_norm in _normalise_strain(c)]
        if len(hits) == 1:
            print(f"Resolved strain query '{query}' -> genome '{hits[0]}' (normalised match)")
            return hits[0]
        if len(hits) > 1:
            raise ValueError(
                f"strain query '{query}' matched multiple genomes: {hits[:10]}"
                + (" ..." if len(hits) > 10 else "")
                + "\n  Re-run with an exact --genome id from this list."
            )

    raise ValueError(
        f"could not resolve '{query}' to any genome. "
        f"Checked exact id, GCF/GCA accession core, and strain name against "
        f"{len(gcols)} genome columns. Example columns: {gcols[:3]}"
    )


def build_presence(
    pa: pd.DataFrame,
    clades: pd.DataFrame,
    family: str,
    clade_id_col: str,
    clade_col: str,
    metadata_cols=None,
    keep_unassigned: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (output_df of genomes WITH the gene, clade_summary over all clades).

    Genomes present in the p/a file but absent from the clade table have no
    clade assignment. By default (keep_unassigned=False) they are dropped from
    the analysis entirely -- the clade table is treated as the definition of
    which genomes are in scope. Set keep_unassigned=True to retain them as an
    'unassigned' group in both the outputs and the figure.
    """
    row = pa.loc[pa["Gene"] == family]
    if row.empty:
        raise ValueError(f"family '{family}' not found in 'Gene' column.")
    row = row.iloc[0]

    gcols = genome_columns(pa, metadata_cols)
    records = []
    for g in gcols:
        ids = cell_ids(row.get(g))
        records.append({"genome": g, "present": int(bool(ids)),
                        "gene_id": ";".join(ids) if ids else pd.NA})
    presence = pd.DataFrame(records)

    clade_map = clades[[clade_id_col, clade_col]].rename(
        columns={clade_id_col: "genome", clade_col: "clade"}
    )
    presence = presence.merge(clade_map, on="genome", how="left")

    missing_mask = presence["clade"].isna()
    n_missing = int(missing_mask.sum())
    if n_missing:
        unassigned_ids = presence.loc[missing_mask, "genome"].tolist()
        if keep_unassigned:
            presence.loc[missing_mask, "clade"] = "unassigned"
            print(f"[warn] {n_missing} genome(s) not in the clade table, kept as "
                  f"'unassigned': {unassigned_ids}", file=sys.stderr)
        else:
            presence = presence.loc[~missing_mask].reset_index(drop=True)
            print(f"[warn] {n_missing} genome(s) not in the clade table, dropped from "
                  f"analysis: {unassigned_ids}\n       (use --keep-unassigned to retain "
                  f"them as an 'unassigned' group)", file=sys.stderr)

    clade_summary = (
        presence.groupby("clade", dropna=False)
        .agg(n_genomes=("genome", "count"), n_with_gene=("present", "sum"))
        .reset_index()
    )
    clade_summary["fraction_with_gene"] = (
        clade_summary["n_with_gene"] / clade_summary["n_genomes"]
    )

    # output_df = only genomes that carry the gene
    output_df = (
        presence.loc[presence["present"] == 1, ["genome", "clade", "gene_id"]]
        .reset_index(drop=True)
    )
    return output_df, clade_summary


def plot_clade_summary(clade_summary: pd.DataFrame, family: str, plot_kind: str = "fraction") -> plt.Figure:
    cs = clade_summary.copy()
    # coerce any NaN clade label to a visible string so it can sit on a
    # categorical axis (fillna before astype -- astype(str) alone can leave
    # NaN as a float in some pandas versions)
    cs["clade"] = cs["clade"].fillna("unassigned").astype(str)
    n = len(cs)
    fig, ax = plt.subplots(figsize=(max(6, 1.1 * n), 5.5))

    n_gene = cs["n_with_gene"].astype(int)
    n_tot = cs["n_genomes"].astype(int)
    if plot_kind == "fraction":
        heights = cs["fraction_with_gene"] * 100
        ylabel, top, off = "Genomes with gene (%)", 100, 2.5
    elif plot_kind == "count":
        heights = n_gene
        ylabel, top = "Number of genomes with gene", max(1, n_gene.max())
        off = 0.03 * top
    else:
        raise ValueError("plot_kind must be 'fraction' or 'count'")

    bars = ax.bar(cs["clade"], heights, color="#4C72B0", edgecolor="white", linewidth=0.5, width=0.65)
    ax.set_ylim(0, top * 1.18)
    for bar, ng, nt in zip(bars, n_gene, n_tot):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + off,
                f"{ng}/{nt}", ha="center", va="bottom", fontsize=9, color="#333")

    ax.set_xlabel("Clade", labelpad=10)
    ax.set_ylabel(ylabel)
    ax.set_title(f"Presence of '{family}' by clade", pad=16, fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, linestyle="--", alpha=0.3)
    ax.set_axisbelow(True)
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
    fig.tight_layout()
    return fig


def run(args: argparse.Namespace) -> None:
    def status(msg: str) -> None:
        if not args.quiet:
            print(msg, flush=True)

    status(f"Loading presence/absence file: {args.pa} ...")
    pa = pd.read_csv(args.pa, dtype=str)
    status(f"File loaded! {pa.shape[0]} gene families x {len(genome_columns(pa))} genomes.")

    # --- inspection mode: just show what a genome's gene ids look like ---
    if args.list_ids:
        if not (args.genome or args.query):
            sys.exit("--list-ids needs a genome: give --genome (exact) or --query (lookup).")
        if args.genome:
            genome = args.genome
        else:
            genome = resolve_genome(pa, args.query)
            status(f"query matched: {genome}")
        print_id_examples(pa, genome, n=args.n_ids)
        return

    missing = [f for f, v in
               [("--clades", args.clades), ("--panaroo-id-col", args.panaroo_id_col),
                ("--clade-col", args.clade_col)] if not v]
    if missing:
        sys.exit(f"missing required argument(s) for analysis: {', '.join(missing)}")
    status(f"Loading clade assignment file: {args.clades} ...")
    clades = pd.read_csv(args.clades)
    status(f"File loaded! {clades.shape[0]} rows.")

    # resolve the target family
    if args.gene_family:
        family = args.gene_family
    else:
        if not args.gene_id:
            sys.exit("Provide either --gene-family, or --gene-id with --genome/--query.")
        if args.genome and args.query:
            sys.exit("Provide only one of --genome (exact) or --query (lookup).")
        if not (args.genome or args.query):
            sys.exit("--gene-id needs a genome: give --genome (exact) or --query (lookup).")

        if args.genome:
            genome = args.genome
        else:
            status(f"Looking up genome for query '{args.query}' ...")
            genome = resolve_genome(pa, args.query)
            status(f"query matched: {genome}")
        check_id_format(pa, genome, args.gene_id)
        family = resolve_family(pa, genome, args.gene_id)
        status(f"Resolved '{args.gene_id}' in '{genome}' -> gene family '{family}'")

    status(f"Mapping clade distribution for gene family '{family}' ...")
    output_df, clade_summary = build_presence(
        pa, clades, family, args.panaroo_id_col, args.clade_col,
        keep_unassigned=args.keep_unassigned,
    )

    status(f"\ngene family '{family}': {output_df.shape[0]} genomes carry the gene")
    if not args.quiet:
        print(clade_summary.to_string(index=False))

    # clades to drop from the FIGURE only (still kept in output_df + printed summary)
    exclude = set(args.exclude_clades) if args.exclude_clades is not None else set(DEFAULT_PLOT_EXCLUDE)
    plot_summary = clade_summary.copy()
    if exclude:
        plot_summary = plot_summary[~plot_summary["clade"].isin(exclude)]
    dropped = sorted(set(clade_summary["clade"].dropna()) - set(plot_summary["clade"].dropna()))
    if dropped:
        status(f"Excluded from figure (still in outputs): {dropped}")

    status("\nPlotting presence per clade ...")
    fig = plot_clade_summary(plot_summary, family, plot_kind=args.plot_kind)
    fig.savefig(args.out_fig, dpi=300, bbox_inches="tight")
    status(f".png file saved to {args.out_fig}")

    if args.out_csv:
        output_df.to_csv(args.out_csv, index=False)
        status(f".csv file saved to {args.out_csv}")

    status("\nDone!")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pa", required=True, help="Panaroo gene_presence_absence.csv")
    p.add_argument("--clades", help="CSV with genome id + clade columns (required unless --list-ids)")
    p.add_argument("--panaroo-id-col", help="column in --clades holding the panaroo genome ids (must match the p/a genome column names)")
    p.add_argument("--clade-col", help="column in --clades holding clade labels")

    g = p.add_argument_group("gene selection (use --gene-family OR --gene-id with a genome)")
    g.add_argument("--gene-family", help="Panaroo gene family (value in the 'Gene' column)")
    g.add_argument("--gene-id", help="a specific gene id within a genome (a leading 'cds-' etc. is matched tolerantly)")
    g.add_argument("--genome", help="EXACT panaroo genome id (a column in the p/a file / panaroo_query_id)")
    g.add_argument("--query", help="loose genome lookup: a GCF/GCA accession or strain name, resolved to a panaroo genome id")

    i = p.add_argument_group("inspection")
    i.add_argument("--list-ids", action="store_true",
                   help="print a sample of real gene ids for the given --genome/--query and exit "
                        "(skips refound/pseudo/len_ and blanks)")
    i.add_argument("--n-ids", type=int, default=20, help="how many example ids to print with --list-ids (default 20)")

    p.add_argument("--out-fig", default="gene_presence_by_clade.png", help="output figure path")
    p.add_argument("--out-csv", default=None, help="optional output CSV of genomes carrying the gene")
    p.add_argument("--plot-kind", choices=["fraction", "count"], default="fraction",
                   help="bar height: fraction of clade (default) or raw count")
    p.add_argument("--exclude-clades", nargs="*", default=None, metavar="CLADE",
                   help="clade label(s) to omit from the FIGURE only (kept in CSV + printed "
                        f"summary). Default drops: {DEFAULT_PLOT_EXCLUDE}. Pass with no values "
                        "to exclude nothing.")
    p.add_argument("--keep-unassigned", action="store_true",
                   help="keep genomes that are in the p/a file but absent from the clade table, "
                        "as an 'unassigned' group in the outputs and figure. Default: drop them "
                        "from the analysis (the clade table defines which genomes are in scope).")
    p.add_argument("--quiet", action="store_true", help="suppress status/progress messages")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
