# Project Overview

## Purpose

This repository contains a single Python tool for processing PALS sMBC
Typhoon `.img` image-plate scans. It locates image files by shot number,
detects one or two 4x4 sMBC arrays, extracts 50 px x 50 px channel
intensities, writes results to Excel, and plots shot-to-shot comparisons.

The main implementation is `smbc_image_tool.py`. User-facing usage examples
live in `README.md`; agent-specific maintenance notes live in `AGENTS.md`.

## Runtime and Dependencies

- Python: tested with Python 3.12.
- Dependencies: `numpy`, `pandas`, `openpyxl`, `matplotlib`, `opencv-python`,
  and `scipy`.
- GUI dependencies:
  - The processing center-adjustment window uses Matplotlib interactivity.
  - The plot-selection GUI uses built-in `tkinter`, imported lazily only for
    `gui-plot`.

Install from the existing requirements file:

```bash
python3 -m pip install -r requirements.txt
```

The project also has minimal `pyproject.toml` metadata and exposes an optional
console entrypoint named `smbc-image-tool` when installed as a package.

## Command Workflows

### `process`

```bash
python3 smbc_image_tool.py process 64433
```

This workflow:

1. Finds the matching `.img` file under the default data root unless `--image`
   is supplied.
2. Reads dimensions from a matching `.inf` file when available, otherwise
   infers square dimensions from file size.
3. Detects byte order automatically unless `--byte-order` is set.
4. Detects channel centers for `layer1`; for shot numbers `64461` and above,
   detects both `layer1` and `layer2`.
5. Opens the manual center-adjustment GUI unless `--no-gui` is set.
6. Extracts ROI statistics and writes Excel sheets.
7. Saves a preview PNG showing centers and ROI boxes.

Important options include `--roi-size`, `--centers-json`,
`--save-centers-json`, `--preview`, `--cmap`, and `--reverse-columns`.

### `plot`

```bash
python3 smbc_image_tool.py plot 64460 64461 64465 \
  --group 64460,64461,64465 \
  --no-show
```

This workflow reads the Excel `intensity_long` sheet, filters by layer, and
saves a comparison PNG. Shot groups share colors and use different markers.
Use `--layer2`, `--range START END`, `--log-y`, `--output`, and `--no-show` as
needed.

### `gui-plot`

```bash
python3 smbc_image_tool.py gui-plot
```

This opens a `tkinter` GUI for the plotting workflow. It reads available shots
from `intensity_long`, lets the user add selected shots into plot groups, then
calls the same `plot_shots()` function used by the CLI. Only shots added to
groups are plotted.

## Data Flow

Default input root:

```text
/Users/zhaoxu/Library/CloudStorage/GoogleDrive-xu.zhao@york.ac.uk/My Drive/Exps/2026 PALS LPI/IP scan
```

Default output Excel:

```text
/Users/zhaoxu/Library/CloudStorage/GoogleDrive-xu.zhao@york.ac.uk/My Drive/Exps/2026 PALS LPI/IP scan/post-processing/sMBC_intensity.xlsx
```

The tool is designed to keep generated data under the experiment
`post-processing` folder by default. Local test outputs can be written to
temporary paths with `--output` or `--preview`.

## Excel Contract

The generated Excel workbook contains:

- `intensity_long`: one row per shot, layer, and channel. Includes channel
  center coordinates, ROI bounds, intensity mean/std/min/max, image metadata,
  ROI size, byte order, and processing timestamp.
- `intensity_wide`: one row per shot. `layer1` channels are `ch01` through
  `ch16`; `layer2` channels are `ch01_2` through `ch16_2`.
- `metadata`: one row per shot with image path, `.inf` path, dimensions, byte
  order, ROI size, processed layers, and timestamp.

Reprocessing the same shot replaces that shot's rows in all generated sheets.

## Implementation Map

- Argument parsing: `parse_args()`.
- Image loading and byte-order detection: `load_img()`, `smoothness_score()`.
- Plate and center detection: `detect_plate_box()`, `detect_centers()`,
  `detect_dual_layer_centers()`, `detect_layer_centers()`.
- Manual center adjustment: `CenterEditor`.
- Measurement extraction and Excel writing: `extract_layer_intensities()`,
  `write_excel()`.
- Preview output: `save_preview()`.
- Plotting: `resolve_plot_shots()`, `resolve_plot_groups()`, `plot_shots()`.
- Plot GUI: `PlotGui`, `run_plot_gui()`.

## Maintenance Notes

- Preserve existing CLI behavior when adding features; many workflows rely on
  command-line compatibility.
- Keep generated files out of Git, especially Excel workbooks, previews, plots,
  caches, and experiment data.
- Keep `tkinter` imports inside the GUI path so headless CLI usage continues to
  work without importing GUI modules.
- Prefer adding focused helper functions over splitting the single-file tool
  prematurely; the current project is intentionally small.
- Run syntax and CLI help checks after changes.
