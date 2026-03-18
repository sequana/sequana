# Copilot Instructions for Sequana

## Project Overview

**Sequana** is a Python library dedicated to bioinformatics, focused on Next-Generation Sequencing (NGS) analysis. It provides standalone tools and Snakemake-based pipelines for quality control, variant calling, genome coverage analysis, taxonomy classification, and transcriptomics.

- **Homepage**: https://sequana.readthedocs.io
- **Repository**: https://github.com/sequana/sequana
- **Python versions**: 3.10, 3.11, 3.12, 3.13

---

## Repository Structure

```
sequana/          # Main Python package
test/             # Tests (pytest)
doc/              # Sphinx documentation
.github/          # GitHub Actions workflows and templates
pyproject.toml    # Project metadata and dependencies (Poetry)
```

Key sub-packages under `sequana/`:
- `enrichment/`     – Gene set enrichment analysis tools
- `kraken/`         – Kraken taxonomy classification tools
- `modules_report/` – HTML report modules
- `multiqc/`        – MultiQC plugin modules
- `plots/`          – Plotting utilities
- `scripts/`        – CLI entry points
- `utils/`          – Utility helpers
- `viz/`            – Visualization tools

---

## Development Setup

Install development dependencies with [Poetry](https://python-poetry.org/):

```bash
pip install poetry
poetry install --with dev,pipelines
```

---

## Coding Conventions

### Code Style
- Formatter: **black** with `--line-length=120`
- Linter: **flake8** (ignore E203, E501, W503, E722; max line length 120)
- Import order: **isort** with `--profile black`
- Pre-commit hooks are configured in `.pre-commit-config.yaml`

Run formatters and linters via pre-commit:
```bash
pre-commit run --all-files
```

### Lazy Imports
Sequana uses a lazy import mechanism to keep import times fast. Common lazy modules are declared in `sequana/lazy.py`. Use `from sequana.lazy import X` for modules like `numpy`, `pandas`, `pysam`, `pylab`, `plotly_express`, `colormap`, `bioservices`, and `gseapy`.

```python
# Prefer this pattern for optional/heavy dependencies:
from sequana.lazy import numpy as np
```

The `LazyImport` class is defined in `sequana/lazyimports.py`.

### General Python Guidelines
- Use Python type hints where appropriate.
- Follow PEP 8 naming conventions (snake_case for functions/variables, CamelCase for classes).
- Keep docstrings concise; use NumPy-style docstrings for public API functions.

---

## Testing

Tests live in the `test/` directory and use **pytest**.

Run the full test suite:
```bash
poetry run pytest --cov-report term-missing --cov=sequana
```

Run a single test file:
```bash
poetry run pytest test/test_<module>.py
```

- Each module `sequana/<name>.py` typically has a corresponding `test/test_<name>.py`.
- Test data is stored in `test/data/`.
- Use `pytest.mark.skipif` or the `skiptravis` fixture from `test/conftest.py` for environment-specific skips.

---

## Pull Request Guidelines

1. Ensure **unit tests**, **documentation**, and **code style** are in order.
2. New functionality should include tests in `test/test_<module>.py`.
3. Contributed code must be compatible with Sequana's BSD license.
4. CI runs on pushes and PRs to `main` and `dev` branches (see `.github/workflows/main.yml`).

---

## Common Pitfalls

- Heavy dependencies (e.g., `pysam`, `bioservices`) must be imported lazily via `sequana.lazy` to avoid slow startup times.
- The package supports Python ≥ 3.10. Avoid using syntax or APIs unavailable in 3.10.
- Snakemake version must be `< 8` (see `pyproject.toml`).
