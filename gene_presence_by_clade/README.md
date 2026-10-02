## Gene presence by clade

`gene_presence_by_clade.py` reports how a chosen gene is distributed across
phylogenetic clades. Given a gene — specified as a gene family, or as a specific
gene ID within a genome — it produces:

- a **bar chart** of the number of genomes per clade that carry the gene, and
- *(optional)* a **CSV** listing the genomes that carry it, with their clade and
  the specific gene ID present.

### Requirements

- Python 3.8+ with `pandas` and `matplotlib`

### Inputs

- **`--pa`** — a Panaroo `gene_presence_absence.csv` (the `.csv`, which contains
  gene IDs per cell — **not** the `.Rtab`).
- **`--clades`** — a metadata table (e.g. `All_Pae_metadata.csv`) with a column
  of genome IDs matching the presence/absence columns and a column of clade
  labels.
- **`--panaroo-id-col`** — the column in `--clades` holding the genome IDs
  (`panaroo_query_id`).
- **`--clade-col`** — the column in `--clades` holding the clade labels.

### Specifying the gene

Use **one** of:

| Flag | Meaning |
|------|---------|
| `--gene-family` | A Panaroo gene family name (from the `Gene` column of the presence/absence file) |
| `--gene-id` + `--genome` | A specific gene ID in an **exact** genome (`panaroo_query_id`) |
| `--gene-id` + `--query` | A specific gene ID in a genome identified **loosely** (accession or strain name) |

### Specifying the genome (for `--gene-id`)

- **`--genome`** — the **exact** `panaroo_query_id` of the genome (found in the
  metadata table).
- **`--query`** — a looser identifier when the exact ID isn't known: a RefSeq
  (`GCF_...`) or GenBank (`GCA_...`) accession, or a strain name (e.g. `PAO1`).
  The script resolves it to a genome automatically.

### Outputs

- **`--out-fig`** — path for the bar chart (PNG).
- **`--out-csv`** *(optional)* — path for the table of genomes carrying the gene,
  with clade and gene ID.

---

### Usage examples

#### 1. A specific gene in a specific genome (exact ID)

`--genome` must be the exact `panaroo_query_id` of the genome (found in
`All_Pae_metadata.csv`).

```bash
python gene_presence_by_clade.py \
  --pa /PATH/TO/gene_presence_absence.csv \
  --clades All_Pae_metadata.csv --panaroo-id-col panaroo_query_id --clade-col clade \
  --genome Pseudomonas_aeruginosa_PAO1_107_GCF_000006765.1 --gene-id NP_249821.1 \
  --out-fig rhlC_distribution_by_clade.png --out-csv rhlC_genomes_and_ids.csv
```

#### 2. A specific gene with partial genome info

If the exact ID isn't known, use `--query` instead of `--genome`, followed by an
accession (`GCF_000006765.1` / `GCA_000006765.1`) or a strain name (`PAO1`); the
script finds the match automatically.

```bash
python gene_presence_by_clade.py \
  --pa /PATH/TO/gene_presence_absence.csv \
  --clades All_Pae_metadata.csv --panaroo-id-col panaroo_query_id --clade-col clade \
  --query GCF_000006765.1 --gene-id NP_249821.1 \
  --out-fig rhlC_distribution_by_clade.png
```

#### 3. Finding the right gene ID

Gene IDs come from the `ID` field of GFF3 files. If unsure which ID to provide,
download the GFF for the genome's accession from NCBI (mind GCF vs GCA). Some IDs
are prefixed with `cds-`; the script accepts the ID with or without it (e.g.
`cds-NP_249821.1` and `NP_249821.1` both work).

Alternatively, use `--list-ids` to print example gene IDs for a genome:

```bash
python gene_presence_by_clade.py \
  --pa /PATH/TO/popPUNK_representatives_gene_presence_absence.csv \
  --clades All_Pae_metadata.csv --panaroo-id-col panaroo_query_id --clade-col clade \
  --query "pseudomonas aeruginosa PAO1" --list-ids
```

#### 4. A specific gene family

If the Panaroo gene family name is known, search by it directly (family names are
in the `Gene` column of the presence/absence files).

```bash
python gene_presence_by_clade.py \
  --pa /PATH/TO/gene_presence_absence.csv \
  --clades All_Pae_metadata.csv --panaroo-id-col panaroo_query_id --clade-col clade \
  --gene-family group_61995 \
  --out-fig exsC_distribution_by_clade.png --out-csv exsC_genomes_and_ids.csv
```
