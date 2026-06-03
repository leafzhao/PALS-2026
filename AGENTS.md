# Agent Notes

## Project Shape

This is a small single-module Python project. The main source file is
`smbc_image_tool.py`; user documentation is in `README.md`; high-level project
context is in `PROJECT_OVERVIEW.md`.

## Working Rules

- Preserve the existing command-line interface unless the user explicitly asks
  for a breaking change.
- Keep `tkinter` imported lazily inside the `gui-plot` path. Normal `process`
  and `plot` commands should not need to import it.
- Do not commit generated experiment data, local Excel workbooks, preview PNGs,
  comparison PNGs, Python caches, or temporary center JSON files.
- Use structured libraries already present in the project (`pandas`,
  `openpyxl`, `numpy`, `opencv-python`, `matplotlib`, `scipy`) instead of ad
  hoc parsing for supported data formats.
- Keep edits scoped. This repository is currently a practical workstation tool,
  not a large package.

## Verification

Run these checks after meaningful code changes:

```bash
python3 -m py_compile smbc_image_tool.py
python3 smbc_image_tool.py --help
python3 smbc_image_tool.py process --help
python3 smbc_image_tool.py plot --help
python3 smbc_image_tool.py gui-plot --help
```

For plotting changes, also run at least one non-interactive plot command with
`--no-show` and an output path outside the repository, for example:

```bash
python3 smbc_image_tool.py plot 64460 64461 \
  --no-show \
  --output /private/tmp/smbc_plot_check.png
```

For metadata changes, validate `pyproject.toml`:

```bash
python3 -c "import pathlib, tomllib; tomllib.loads(pathlib.Path('pyproject.toml').read_text())"
```

## Git Hygiene

- Commit source, docs, dependency files, and project metadata.
- Leave `__pycache__/` and generated outputs untracked.
- If the working tree is dirty, inspect changes before committing and avoid
  reverting user work.
