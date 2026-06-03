# PALS sMBC Image Tool

This tool processes Typhoon `.img` files for the sMBC diagnostic array. It
accepts a shot number, locates the matching image, detects sMBC 4x4 regions
while excluding right-side long rectangular strips, lets the user manually
adjust centers, averages a 50 px x 50 px ROI for each channel, writes the
result to Excel, and plots shot-to-shot comparisons. Starting at shot `64461`,
the tool processes two 4x4 layers: the left array is `layer1`, and the right
array is `layer2`.

## Project Docs

- `PROJECT_OVERVIEW.md`: architecture, data flow, Excel sheet contract, and
  maintenance notes.
- `AGENTS.md`: repo-specific guidance for future coding-agent work.
- `pyproject.toml`: minimal Python project metadata and optional console
  script entrypoint.

## Install

```bash
python3 -m pip install -r requirements.txt
```

The current workstation already has the required packages installed.

## Process One Shot

```bash
python3 smbc_image_tool.py process 64433
```

Default input root:

```text
/Users/zhaoxu/Library/CloudStorage/GoogleDrive-xu.zhao@york.ac.uk/My Drive/Exps/2026 PALS LPI/IP scan
```

Default output:

```text
/Users/zhaoxu/Library/CloudStorage/GoogleDrive-xu.zhao@york.ac.uk/My Drive/Exps/2026 PALS LPI/IP scan/post-processing/sMBC_intensity.xlsx
```

All default processed outputs are written under the `post-processing` folder,
which is created automatically when needed.

The GUI shows the detected centers, channel labels, and each 50 px x 50 px ROI.
Drag a channel marker or box to adjust one center. Press `m` to toggle moving
the whole grid. Arrow keys move the selected channel by 1 px; Shift+arrow moves
by 10 px. Press Enter to accept and write the Excel file. Press `q` or Escape to
abort.

For shots `64461` and later, preview labels use `1..16` for `layer1` and
`1_2..16_2` for `layer2`.

For automatic batch mode without manual adjustment:

```bash
python3 smbc_image_tool.py process 64433 --no-gui
```

Useful options:

```bash
python3 smbc_image_tool.py process 64433 \
  --roi-size 50 \
  --excel ./sMBC_intensity.xlsx \
  --preview ./sMBC_64433_preview.png \
  --cmap jet \
  --save-centers-json ./sMBC_64433_centers.json
```

Preview images and the manual adjustment window use the `jet` colormap by
default and include an intensity colorbar. Change it with `--cmap` or
`--preview-cmap`, for example `--cmap gray`, `--cmap viridis`, or `--cmap turbo`.

Channel numbering defaults to bottom-left row-major order: bottom row is
channels 1-4 from left to right, the next row is 5-8, and the top row is 13-16.
Use `--reverse-columns` if a dataset needs right-to-left numbering inside each
row.

## Plot Shot Comparisons

```bash
python3 smbc_image_tool.py plot 64433 64434 64435
```

This reads `intensity_long` from the Excel file and saves a PNG with channel on
the x-axis and 50 px x 50 px mean intensity on the y-axis. The plot defaults to
`layer1`; use `--layer2` or `-layer2` to compare the second layer. The plot
window opens by default after saving. Comparison PNG names include the layer,
for example `sMBC_compare_layer1_64433_64434.png`; if that file already exists,
the tool creates `_001`, then `_002`, and so on.

```bash
python3 smbc_image_tool.py plot 64433 64434 \
  --excel ./sMBC_intensity.xlsx \
  --output ./sMBC_compare.png
```

Compare a closed shot range using shots already present in the Excel file:

```bash
python3 smbc_image_tool.py plot --range 64433 64461
python3 smbc_image_tool.py plot --range 64433 64461 --layer2
```

Shots missing from the Excel file are skipped with a warning. When plotting
`layer2`, shots without `layer2` data are also skipped with a warning.

Use a logarithmic intensity axis with:

```bash
python3 smbc_image_tool.py plot 64433 64434 --log-y
python3 smbc_image_tool.py plot --range 64461 64470 --layer2 --log-y
```

Group shots so each group shares one color:

```bash
python3 smbc_image_tool.py plot 64433 64434 64461 \
  --group 64433,64434 \
  --group 64461
```

Grouped shots share a color and use different markers. Shots not listed in any
`--group` are still plotted as one-shot groups with their own colors.

For batch plotting without opening a window:

```bash
python3 smbc_image_tool.py plot 64433 64434 --no-show
```

### Plot with the GUI

Use the GUI when a comparison needs many shots or color groups:

```bash
python3 smbc_image_tool.py gui-plot
```

The GUI loads shots from the `intensity_long` Excel sheet. Select one or more
shots in the left list and click `Add group >`; each group shares one plot color
and uses different markers, matching the `--group` command-line behavior. Only
shots added to groups are plotted. Use the layer, log-y, show-plot, and optional
output path controls, then click `Run plot`.

You can start the GUI with a different Excel file:

```bash
python3 smbc_image_tool.py gui-plot --excel ./sMBC_intensity.xlsx
```

## Excel Sheets

- `intensity_long`: one row per shot, layer, and channel, including center
  coordinates, ROI bounds, mean, standard deviation, min, and max.
- `intensity_wide`: one row per shot with `layer1` columns `ch01` to `ch16`;
  `layer2` columns use `_2`, such as `ch01_2` to `ch16_2`.
- `metadata`: image path, dimensions, byte order, ROI size, processed layers,
  and processing time.

Processing the same shot again overwrites that shot's rows in all generated
sheets.
