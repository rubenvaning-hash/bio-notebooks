# bio-notebooks

Data analysis notebooks in Python, built alongside the Bioinformatics master at UGent.

The aim is applied rather than theoretical: load real biological data, clean it,
plot it, and say something defensible about what it shows.

## Setup

Dependencies are managed with [uv](https://github.com/astral-sh/uv). To rebuild
the exact environment on any machine:

```
uv sync
```

This reads `uv.lock` and creates a local `.venv` with the pinned versions,
including the Python interpreter itself. Nothing is installed system-wide.

## Notebooks

| Notebook | Dataset | What it does |
|---|---|---|
| `start.ipynb` | `seaborn` penguins (344 rows) | Environment check. Shape, first rows, and species means. |

## Notes

Notebooks are committed **with their outputs** on purpose, so results are
readable on GitHub without cloning and running anything.
