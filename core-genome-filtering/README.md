## Core-genome filtering for alignments

A pipeline to take a core-genome generated from Panaroo output, remove gene families dominated by paralog and pseudogene flags, clean up "spurious" flags, and generate nucleotide FASTAs per gene family that can be used as input for alignment tools like MAFFT. 


### Overview

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `plot_flagged_distribution.py` | Plot paralog / pseudo / refound flag distributions to choose thresholds |
| 2 | `filter_core_genome.py` | Remove flag-dominated families, blank spurious flagged cells, drop low-fill genomes |
| 3 | `write_gene_fastas.py` | Write one nucleotide FASTA per gene family from `gene_data.csv` |
| – | `subset_gene_data.py` | (companion) Write a core-only subset of the huge `gene_data.csv` |
| – | `verify_gene_data_subset.py` | (companion) Validate a `gene_data.csv` subset for completeness |

### Requirements

- Python 3.8+
- `pandas`, `numpy`, `matplotlib`

```bash
pip install pandas numpy matplotlib
```

### Inputs

- **`core4079_gene_presence_absence.csv`** — the core pan-genome presence/absence
  table. Must be the **`.csv`** (contains gene IDs per cell), **not** the `.Rtab`
  (0/1 only). First 3 columns are metadata (`Gene`, `Non-unique Gene name`,
  `Annotation`); remaining columns are genomes.
- **`core4079_gene_data.csv`** — the Panaroo `merge` output mapping every gene ID
  (`annotation_id`) to its nucleotide sequence (`dna_sequence`).
- *(optional)* a metadata table (e.g. `All_Pae_metadata.csv`) and a subset list
  of genomes, if restricting the analysis to a subset.

> **Recommended starting point:** begin from the provided
> `core4079_gene_presence_absence.csv` (the 98%-presence core).
>
> For a **de novo** definition of the core genome, or to **subset from the full
> pan-genome**, download the full pan-genome files
> (`gene_presence_absence.csv`, `gene_presence_absence.Rtab`, `gene_data.csv`),
> define/filter your desired core from those first, then use the result as input
> to Step 1.

---

### Step 1 — Inspect flag distributions and choose thresholds

Start from the 98%-presence core (`core4079_gene_presence_absence.csv`) and plot
the distribution of gene families containing paralogs and pseudo/refound flags,
along with the number of genomes carrying each flag.

```bash
python plot_flagged_distribution.py \
  --pa /PATH/TO/core4079_gene_presence_absence.csv \
  --out-dir /PATH/TO/OUTDIR
```

**Analyzing a subset of genomes:** add `--subset`, a text file with one genome
ID per line. IDs may be RefSeq (`GCF_...`) or GenBank (`GCA_...`) accessions, or
exact `panaroo_query_id` values; each line can use any recognizable format. They
are resolved against the presence/absence column names.

```bash
python plot_flagged_distribution.py \
  --pa /PATH/TO/core4079_gene_presence_absence.csv \
  --subset /PATH/TO/popPUNK_cluster_representatives.txt \
  --out-dir /PATH/TO/OUTDIR
```

**Outputs:** PNG files of the paralog and pseudo/refound flag distributions.
(`refound` = genes not annotated in the original GFF but recovered by Panaroo's
BLAST step.)

**Finding the right thresholds:** paralogs can confuse phylogenetic
analysis. The paralog plot shows a large peak at low genome counts (many gene
families with paralogs in only a handful of genomes) and occasional bars far out
at hundreds to thousands of genomes (1–2 families each). The far-right bars represent
**genuine** gene duplications; the low peak is mostly **spurious** —
fragmentation or upstream-pipeline artifacts. Step 2 removes the genuine
multi-copy families entirely, and blanks the individual spurious calls.

The script prints, for each flag type, how many families a given threshold would
drop — use that table plus the plots to pick your Step 2 thresholds.

---

### Step 2 — Filter the core genome


1. **Remove gene families** where paralog / pseudo flags appear in a majority of
   genomes (genuine paralogous or pseudogene families).
2. **Blank remaining spurious flagged cells** — set those specific cells to
   empty. Downstream steps treat a blank cell as *missing information*, not as an
   absent gene.
3. **Drop low-content genomes** — remove genomes carrying less than a chosen
   fraction of the core genes after cleaning.

```bash
python filter_core_genome.py \
  --pa /PATH/TO/core4079_gene_presence_absence.csv \
  --paralog-threshold 10 --pseudo-threshold 1 --min-fill 80 \
  --out /PATH/TO/core_gene_presence_absence_filtered.csv
```

- Thresholds are **percentages** (`2` → 2% of genomes).
- `--min-fill` (default `80`) sets the genome-completeness cutoff.
- `refound` flags are left untouched.

**Outputs** (written next to `--out`):

- the filtered presence/absence CSV,
- `*.filter_log.txt` — thresholds used and the names of every dropped gene family and
  genome (for reproducibility),
- `*.fill_per_genome.png` — a before/after QC plot of % core genes filled per
  genome.

---

### Step 3 — Write per-gene FASTA files

Build one nucleotide FASTA per gene family, ready for input into alignment tool. Sequences come
directly from `core4079_gene_data.csv` (keyed by the gene IDs in the filtered
presence/absence cells) — **no genome downloads required**.

```bash
python write_gene_fastas.py \
  --pa /PATH/TO/core_gene_presence_absence_filtered.csv \
  --gene-data /PATH/TO/core4079_gene_data.csv \
  --out-dir /PATH/TO/mafft_input_fasta
```

Each FASTA holds one sequence per genome that carries the gene, with the genome
ID as the FASTA header.

> **Tip:** point `--gene-data` at the core-only subset (`core4079_gene_data.csv`,
> from `subset_gene_data.py` below) instead of the full `gene_data.csv` for a
> faster run.

---

### Companion — subset `gene_data.csv` to the core (`core4079_gene_data.csv`)

`gene_data.csv` is a large file. This writes a smaller subset (example for core4079 `core4079_gene_data.csv` here) containing only the rows for the gene families in a given core presence/absence table.

```bash
python subset_gene_data.py \
  --pa /PATH/TO/core4079_gene_presence_absence.csv \
  --gene-data /PATH/TO/gene_data.csv \
  --out /PATH/TO/core4079_gene_data.csv \
  --write-missing /PATH/TO/missing_ids.txt
```

It streams `gene_data.csv` and reports coverage (`X/Y core ids matched a gene_data row`). Any core ID with no matching row is written to `--write-missing`.



#### Verifying a subset

To confirm a subset is complete (e.g. after a long run), compare it against the
core presence/absence it was built from:

```bash
python verify_gene_data_subset.py \
  --pa /PATH/TO/core4079_gene_presence_absence.csv \
  --subset /PATH/TO/core4079_gene_data.csv \
  --progress
```

It reports coverage (missing IDs), empty-sequence rows, and malformed/truncated
rows, and exits non-zero if any hard check fails. Given the note above, expect
the only "missing" IDs to be `*_pseudo`.

---
