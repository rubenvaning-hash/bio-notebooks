# bio-notebooks

Data analysis notebooks in Python, built alongside the Bioinformatics master at UGent.

The aim is applied rather than theoretical: load real biological data, clean it,
plot it, and say something defensible about what it shows.

## Featured: dexamethasone in airway smooth muscle cells

**Notebook:** [`04-airway.ipynb`](04-airway.ipynb)

**Question.** Which genes change when airway smooth muscle cells are treated with
dexamethasone, a glucocorticoid used against asthma?

**Data.** RNA-seq counts from Himes et al. (2014), *PLoS ONE* 9(6): e99625,
GEO accession [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778).
Four cell lines, each measured once untreated and once treated: 38,694 genes x 8 samples.

**Answer.** With only 4 pairs, a per-gene paired t-test finds just 1 significant gene
after Benjamini-Hochberg correction, even though the paper's headline gene, CRISPLD2,
goes up about 5x (p = 0.004). Borrowing spread information across all genes
(a moderated t-statistic) moves CRISPLD2 and three other known glucocorticoid-response
genes sharply up the ranking.

**How it gets there:**

1. **Filter.** 13,436 genes are zero in every sample. Removing only those is not enough:
   11 of the first 15 "hits" were genes with about one read, where all four differences
   were identical, so the spread was 0 and the t-test returned p = 0. Final filter:
   mean count >= 10, leaving 15,285 genes.
2. **Test.** log2(count + 1), then a paired t-test per gene (each cell line is its own
   control), then Benjamini-Hochberg written by hand and checked against scipy.
3. **Result.** 1 gene survives at q < 0.05. A volcano plot shows why: the BH cutoff sits
   at -log10(p) = 5.65, CRISPLD2 at 2.4.
4. **Moderated t.** Adding s0 (the median spread across all genes) to every gene's spread
   stops genes with a suspiciously tiny spread from dominating the list.

| Gene | Rank, normal t | Rank, moderated t |
|---|---|---|
| CRISPLD2 | 156 | 82 |
| FKBP5 | 68 | 20 |
| KLF15 | 75 | 18 |
| DUSP1 | 95 | 40 |

**Why it matters.** This is the problem that limma and DESeq2 solve: with few samples,
a gene's own spread is too unreliable to test it alone.

**Limits.** n = 4 pairs. The moderated t gives ranks, not p-values. The counts are loaded
from a course copy of the dataset (bioboot BIMM143), not processed from the raw GEO files.

## Notebooks

| Notebook | Dataset | Finding |
|---|---|---|
| `start.ipynb` | `seaborn` penguins (344 rows) | Environment check: shape, first rows, species means. |
| [`01-first-steps.ipynb`](01-first-steps.ipynb) | penguins | Boolean filters silently drop rows with missing values. |
| [`02-one-figure.ipynb`](02-one-figure.ipynb) | penguins | Bill length alone cannot tell Chinstrap from Gentoo; adding bill depth separates all three species. |
| [`03-multiple-testing.ipynb`](03-multiple-testing.ipynb) | simulated, 10,000 genes | Without correction, 38% of the "significant" list is false. Bonferroni keeps only 3 real effects; Benjamini-Hochberg (written by hand, identical to scipy) keeps 82 with 0 false. |
| [`04-airway.ipynb`](04-airway.ipynb) | airway, GSE52778 | See above: with n = 4, a per-gene t-test misses the paper's headline gene. |

## Setup

Dependencies are managed with [uv](https://github.com/astral-sh/uv). To rebuild
the exact environment on any machine:

```
uv sync
```

This reads `uv.lock` and creates a local `.venv` with the pinned versions,
including the Python interpreter itself. Nothing is installed system-wide.
Then open a notebook in VS Code and select the `.venv` kernel.
`04-airway.ipynb` downloads its data on first run, so it needs an internet connection.

## Notes

Notebooks are committed **with their outputs** on purpose, so results are
readable on GitHub without cloning and running anything.

Notebooks `02` and `04` end with a short note on how AI was used: what was typed
from suggested code and what was written independently. This README was drafted
by Claude from the notebooks' findings.
