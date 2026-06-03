#!/usr/bin/env python3
"""Process sMBC 4x4 image-plate data from Typhoon .img files."""

from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import hashlib
import io
import json
import math
import os
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "pals_matplotlib"))

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.signal import find_peaks


DEFAULT_DATA_ROOT = Path(
    "/Users/zhaoxu/Library/CloudStorage/GoogleDrive-xu.zhao@york.ac.uk/"
    "My Drive/Exps/2026 PALS LPI/IP scan"
)
DEFAULT_OUTPUT_DIR = DEFAULT_DATA_ROOT / "post-processing"
DEFAULT_EXCEL = DEFAULT_OUTPUT_DIR / "sMBC_intensity.xlsx"
DEFAULT_ROI_SIZE = 50
DUAL_LAYER_START_SHOT = 64461
CHANNELS = list(range(1, 17))


@dataclass(frozen=True)
class ImageInfo:
    path: Path
    inf_path: Path | None
    width: int
    height: int
    dtype: str


@dataclass(frozen=True)
class PlateBox:
    left: int
    top: int
    width: int
    height: int

    @property
    def right(self) -> int:
        return self.left + self.width

    @property
    def bottom(self) -> int:
        return self.top + self.height


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract 50x50 mean intensities from sMBC 4x4 .img data."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    process = subparsers.add_parser("process", help="process one shot")
    process.add_argument("shot", help="shot number, e.g. 64433")
    process.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    process.add_argument("--excel", type=Path, default=DEFAULT_EXCEL)
    process.add_argument("--image", type=Path, help="explicit .img file path")
    process.add_argument("--roi-size", type=int, default=DEFAULT_ROI_SIZE)
    process.add_argument(
        "--byte-order",
        choices=("auto", "big", "little"),
        default="auto",
        help="raw .img byte order; Typhoon files here are normally big endian",
    )
    process.add_argument(
        "--no-gui",
        action="store_true",
        help="skip manual adjustment window and use automatic centers",
    )
    process.add_argument(
        "--preview",
        type=Path,
        help="save overlay preview PNG; default is next to the Excel file",
    )
    process.add_argument(
        "--preview-cmap",
        "--cmap",
        default="jet",
        help="Matplotlib colormap for preview and manual GUI; default is jet",
    )
    process.add_argument(
        "--centers-json",
        type=Path,
        help="load center coordinates from a JSON file instead of auto detection",
    )
    process.add_argument(
        "--save-centers-json",
        type=Path,
        help="save final center coordinates to JSON",
    )
    process.add_argument(
        "--reverse-columns",
        action="store_true",
        help="number each row right-to-left instead of left-to-right",
    )

    plot = subparsers.add_parser("plot", help="plot shots from the Excel file")
    plot.add_argument("shots", nargs="*", help="shot numbers to compare")
    plot.add_argument(
        "--range",
        dest="shot_range",
        nargs=2,
        metavar=("START", "END"),
        help="inclusive shot range to compare, using shots already present in Excel",
    )
    plot.add_argument("--excel", type=Path, default=DEFAULT_EXCEL)
    plot.add_argument(
        "--output",
        type=Path,
        help="output PNG path; a suffix is added if the file already exists",
    )
    plot.add_argument(
        "--show",
        action="store_true",
        default=True,
        help="open an interactive plot window after saving; enabled by default",
    )
    plot.add_argument(
        "--no-show",
        action="store_false",
        dest="show",
        help="save the plot without opening a window",
    )
    plot.add_argument(
        "--log-y",
        action="store_true",
        help="use a logarithmic y-axis for intensity",
    )
    plot.add_argument(
        "--group",
        action="append",
        default=[],
        help="comma-separated shots that should share one plot color; repeat for multiple groups",
    )
    layer_group = plot.add_mutually_exclusive_group()
    layer_group.add_argument(
        "--layer1",
        "-layer1",
        action="store_const",
        const="layer1",
        dest="layer",
        default="layer1",
        help="plot layer1 data; this is the default",
    )
    layer_group.add_argument(
        "--layer2",
        "-layer2",
        action="store_const",
        const="layer2",
        dest="layer",
        help="plot layer2 data",
    )

    gui_plot = subparsers.add_parser("gui-plot", help="open a GUI for plotting shot comparisons")
    gui_plot.add_argument(
        "--excel",
        type=Path,
        default=DEFAULT_EXCEL,
        help="initial Excel file to load",
    )

    return parser.parse_args()


def read_inf_dimensions(inf_path: Path) -> tuple[int, int] | None:
    if not inf_path.exists():
        return None
    numeric: list[int] = []
    for raw_line in inf_path.read_text(errors="ignore").splitlines():
        line = raw_line.strip()
        if line.isdigit():
            numeric.append(int(line))
    for idx, value in enumerate(numeric):
        if value == 16 and idx + 2 < len(numeric):
            width, height = numeric[idx + 1], numeric[idx + 2]
            if width > 0 and height > 0:
                return width, height
    return None


def locate_image(shot: str, data_root: Path, explicit_image: Path | None) -> Path:
    if explicit_image is not None:
        if not explicit_image.exists():
            raise FileNotFoundError(f"Image file not found: {explicit_image}")
        return explicit_image

    matches = sorted(data_root.rglob(f"{shot}*.img"))
    if not matches:
        raise FileNotFoundError(f"No .img file found for shot {shot} under {data_root}")
    exact_prefix = [path for path in matches if path.name.startswith(f"{shot}-")]
    matches = exact_prefix or matches
    phosphor = [path for path in matches if "phosphor" in path.name.lower()]
    matches = phosphor or matches
    if len(matches) > 1:
        hashes: dict[str, list[Path]] = {}
        for path in matches:
            hashes.setdefault(file_sha256(path), []).append(path)
        if len(hashes) == 1:
            return matches[0]
        groups = []
        for digest, paths in hashes.items():
            groups.append(digest + ":\n" + "\n".join(f"  {path}" for path in paths))
        options = "\n".join(groups)
        raise RuntimeError(
            f"Multiple different .img files found for shot {shot}; "
            f"use --image to choose one:\n{options}"
        )
    return matches[0]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def infer_square_dimensions(path: Path) -> tuple[int, int]:
    pixels = path.stat().st_size // 2
    side = int(math.sqrt(pixels))
    if side * side != pixels:
        raise ValueError(
            f"Cannot infer dimensions from {path}; provide a matching .inf file"
        )
    return side, side


def smoothness_score(data: np.ndarray) -> float:
    sample = data.astype(np.float32)[::4, ::4]
    diffs = np.mean(np.abs(np.diff(sample, axis=0))) + np.mean(
        np.abs(np.diff(sample, axis=1))
    )
    return float(diffs / (np.std(sample) + 1e-9))


def load_img(path: Path, byte_order: str = "auto") -> tuple[np.ndarray, ImageInfo]:
    inf_path = path.with_suffix(".inf")
    dims = read_inf_dimensions(inf_path)
    width, height = dims if dims is not None else infer_square_dimensions(path)

    dtype_map = {"big": ">u2", "little": "<u2"}
    if byte_order in dtype_map:
        dtype = dtype_map[byte_order]
    else:
        candidates: dict[str, np.ndarray] = {}
        scores: dict[str, float] = {}
        for dtype_name in (">u2", "<u2"):
            arr = np.fromfile(path, dtype=np.dtype(dtype_name)).reshape(height, width)
            candidates[dtype_name] = arr
            scores[dtype_name] = smoothness_score(arr)
        dtype = min(scores, key=scores.get)
        image = candidates[dtype]
        return image, ImageInfo(path, inf_path if inf_path.exists() else None, width, height, dtype)

    image = np.fromfile(path, dtype=np.dtype(dtype)).reshape(height, width)
    return image, ImageInfo(path, inf_path if inf_path.exists() else None, width, height, dtype)


def display_limits(image: np.ndarray) -> tuple[float, float]:
    lo, hi = np.percentile(image, [1, 99.5])
    if hi <= lo:
        hi = float(np.max(image))
        lo = float(np.min(image))
    return float(lo), float(hi)


def make_unique_path(path: Path) -> Path:
    if not path.exists():
        return path
    for index in range(1, 10000):
        candidate = path.with_name(f"{path.stem}_{index:03d}{path.suffix}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"Could not find a free output filename for {path}")


def compare_plot_name(
    shots: list[str],
    layer: str,
    log_y: bool = False,
    grouped: bool = False,
) -> str:
    joined = "_".join(shots)
    safe = "".join(char if char.isalnum() or char in ("-", "_") else "_" for char in joined)
    if len(safe) > 120:
        safe = safe[:120].rstrip("_")
    suffix = ""
    if grouped:
        suffix += "_grouped"
    if log_y:
        suffix += "_logy"
    return f"sMBC_compare_{layer}_{safe}{suffix}.png"


def parse_shot_number(shot: str) -> int | None:
    try:
        return int(str(shot))
    except ValueError:
        return None


def is_dual_layer_shot(shot: str) -> bool:
    shot_number = parse_shot_number(shot)
    return shot_number is not None and shot_number >= DUAL_LAYER_START_SHOT


def channel_label(layer: str, channel: int) -> str:
    if layer == "layer1":
        return str(channel)
    suffix = layer.removeprefix("layer")
    return f"{channel}_{suffix}"


def layer_sort_key(layer: str) -> tuple[int, str]:
    if layer.startswith("layer"):
        try:
            return int(layer.removeprefix("layer")), layer
        except ValueError:
            pass
    return 999, layer


def union_plate_box(boxes: list[PlateBox]) -> PlateBox:
    left = min(box.left for box in boxes)
    top = min(box.top for box in boxes)
    right = max(box.right for box in boxes)
    bottom = max(box.bottom for box in boxes)
    return PlateBox(left, top, right - left, bottom - top)


def box_from_grid(
    xs: list[float],
    ys_bottom_up: list[float],
    image_shape: tuple[int, int],
) -> PlateBox:
    height, width = image_shape
    sorted_xs = sorted(xs)
    sorted_ys = sorted(ys_bottom_up)
    dx = float(np.median(np.diff(sorted_xs))) if len(sorted_xs) > 1 else width * 0.1
    dy = float(np.median(np.diff(sorted_ys))) if len(sorted_ys) > 1 else height * 0.1
    left = max(0, int(round(min(sorted_xs) - 0.35 * dx)))
    right = min(width, int(round(max(sorted_xs) + 0.35 * dx)))
    top = max(0, int(round(min(sorted_ys) - 0.45 * dy)))
    bottom = min(height, int(round(max(sorted_ys) + 0.40 * dy)))
    return PlateBox(left, top, right - left, bottom - top)


def boxes_from_centers(
    centers_by_layer: dict[str, dict[int, tuple[float, float]]],
    image_shape: tuple[int, int],
) -> dict[str, PlateBox]:
    boxes: dict[str, PlateBox] = {}
    for layer, centers in centers_by_layer.items():
        xs = [center[0] for center in centers.values()]
        ys = [center[1] for center in centers.values()]
        boxes[layer] = box_from_grid(xs, ys, image_shape)
    return boxes


def detect_plate_box(image: np.ndarray) -> PlateBox:
    height, width = image.shape
    work = image.astype(np.float32)
    smooth = cv2.GaussianBlur(work, (0, 0), 12)
    bg = float(np.percentile(smooth, 5))
    hi = float(np.percentile(smooth, 95))
    threshold = bg + 0.12 * (hi - bg)
    mask = (smooth > threshold).astype(np.uint8) * 255
    kernel = np.ones((31, 31), np.uint8)
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
    count, _labels, stats, _centroids = cv2.connectedComponentsWithStats(mask)

    square_candidates: list[tuple[int, int, int, int, int]] = []
    candidates: list[tuple[int, int, int, int, int]] = []
    image_area = image.shape[0] * image.shape[1]
    for idx in range(1, count):
        area = int(stats[idx, cv2.CC_STAT_AREA])
        x = int(stats[idx, cv2.CC_STAT_LEFT])
        y = int(stats[idx, cv2.CC_STAT_TOP])
        w = int(stats[idx, cv2.CC_STAT_WIDTH])
        h = int(stats[idx, cv2.CC_STAT_HEIGHT])
        if area < image_area * 0.02:
            continue
        aspect = w / max(h, 1)
        center_x = x + w / 2
        if 0.70 <= aspect <= 1.30 and center_x < width * 0.75:
            square_candidates.append((area, x, y, w, h))
        if 0.45 <= aspect <= 2.2:
            candidates.append((area, x, y, w, h))

    if square_candidates:
        _area, x, y, w, h = max(square_candidates, key=lambda item: item[0])
        return PlateBox(x, y, w, h)

    if not candidates:
        margin_x = int(width * 0.12)
        margin_y = int(height * 0.12)
        return PlateBox(margin_x, margin_y, width - 2 * margin_x, height - 2 * margin_y)

    _area, x, y, w, h = max(candidates, key=lambda item: item[0])
    return PlateBox(x, y, w, h)


def set_plate_view(ax, image_shape: tuple[int, int], plate_box: PlateBox, margin: int = 80) -> None:  # noqa: ANN001
    height, width = image_shape
    right_margin = min(margin, max(20, int(plate_box.width * 0.03)))
    ax.set_xlim(max(0, plate_box.left - margin), min(width, plate_box.right + right_margin))
    ax.set_ylim(min(height, plate_box.bottom + margin), max(0, plate_box.top - margin))


def find_profile_peaks(
    profile: np.ndarray,
    expected: int,
    min_distance: int,
    prominence_scale: float = 0.25,
) -> list[int]:
    smooth = cv2.GaussianBlur(profile.astype(np.float32).reshape(1, -1), (0, 0), 35).ravel()
    prominence = max(float(np.std(smooth)) * prominence_scale, 1.0)
    peaks, properties = find_peaks(smooth, distance=max(min_distance, 1), prominence=prominence)
    if len(peaks) == 0:
        return []
    prominences = properties.get("prominences", np.zeros(len(peaks)))
    ranked = sorted(zip(peaks, prominences), key=lambda item: item[1], reverse=True)
    selected = sorted(int(peak) for peak, _prom in ranked[:expected])
    return selected


def fallback_grid_from_box(box: PlateBox) -> tuple[list[float], list[float]]:
    x_margin = 0.095 * box.width
    y_top_margin = 0.12 * box.height
    y_bottom_margin = 0.10 * box.height
    xs = np.linspace(box.left + x_margin, box.right - x_margin, 4).tolist()
    ys_bottom_up = np.linspace(box.bottom - y_bottom_margin, box.top + y_top_margin, 4).tolist()
    return xs, ys_bottom_up


def extrapolate_rows_from_peaks(peaks_y: list[int], box: PlateBox) -> list[float]:
    _xs, fallback_rows = fallback_grid_from_box(box)
    if not peaks_y:
        return fallback_rows

    peaks = sorted(peaks_y, reverse=True)
    bottom = float(peaks[0])
    if len(peaks) >= 2:
        diffs = [peaks[i] - peaks[i + 1] for i in range(len(peaks) - 1)]
        spacing = float(np.median([diff for diff in diffs if diff > 0] or [box.height / 4]))
    else:
        top_guess = box.top + 0.12 * box.height
        spacing = max((bottom - top_guess) / 3.0, box.height / 5.0)

    rows = [bottom - idx * spacing for idx in range(4)]
    top_limit = box.top + 0.04 * box.height
    bottom_limit = box.bottom - 0.03 * box.height
    if rows[-1] < top_limit or rows[0] > bottom_limit or spacing <= 0:
        return fallback_rows
    return rows


def refine_center(image: np.ndarray, center: tuple[float, float], roi_size: int) -> tuple[float, float]:
    x0, y0 = center
    radius = max(roi_size * 2, 100)
    x1 = max(0, int(round(x0 - radius)))
    x2 = min(image.shape[1], int(round(x0 + radius)))
    y1 = max(0, int(round(y0 - radius)))
    y2 = min(image.shape[0], int(round(y0 + radius)))
    if x2 <= x1 or y2 <= y1:
        return center

    patch = image[y1:y2, x1:x2].astype(np.float32)
    background = cv2.GaussianBlur(patch, (0, 0), max(roi_size, 25))
    residual = patch - background
    median = float(np.median(residual))
    high = float(np.percentile(residual, 92))
    if high <= median:
        return center
    mask = (residual > median + 0.45 * (high - median)).astype(np.uint8)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))
    count, labels, stats, centroids = cv2.connectedComponentsWithStats(mask)
    best_idx = None
    best_score = float("inf")
    local_x = x0 - x1
    local_y = y0 - y1
    for idx in range(1, count):
        area = int(stats[idx, cv2.CC_STAT_AREA])
        if area < roi_size * roi_size * 0.08:
            continue
        cx, cy = centroids[idx]
        distance = math.hypot(float(cx) - local_x, float(cy) - local_y)
        if distance < best_score:
            best_idx = idx
            best_score = distance
    if best_idx is None or best_score > radius * 0.75:
        return center
    cx, cy = centroids[best_idx]
    return float(x1 + cx), float(y1 + cy)


def detect_centers(
    image: np.ndarray,
    roi_size: int,
    reverse_columns: bool = False,
) -> tuple[dict[int, tuple[float, float]], PlateBox]:
    box = detect_plate_box(image)
    roi = image[box.top : box.bottom, box.left : box.right].astype(np.float32)
    xs_fallback, ys_fallback = fallback_grid_from_box(box)

    bottom_start = int(box.height * 0.78)
    bottom_stop = int(box.height * 0.97)
    bottom_band = roi[bottom_start:bottom_stop, :]
    x_peaks = find_profile_peaks(
        bottom_band.mean(axis=0),
        expected=4,
        min_distance=max(int(box.width / 6), roi_size * 2),
    )
    xs = [float(box.left + peak) for peak in x_peaks] if len(x_peaks) == 4 else xs_fallback

    left_stop = max(int(box.width * 0.24), roi_size * 3)
    left_band = roi[:, :left_stop]
    y_peaks = find_profile_peaks(
        left_band.mean(axis=1),
        expected=4,
        min_distance=max(int(box.height / 7), roi_size * 2),
    )
    ys_bottom_up = extrapolate_rows_from_peaks([box.top + peak for peak in y_peaks], box)

    centers: dict[int, tuple[float, float]] = {}
    channel = 1
    row_xs = list(reversed(xs)) if reverse_columns else xs
    for y in ys_bottom_up:
        for x in row_xs:
            centers[channel] = refine_center(image, (x, y), roi_size)
            channel += 1
    return centers, box


def detect_dual_layer_array_box(image: np.ndarray) -> PlateBox:
    height, width = image.shape
    work = image.astype(np.float32)
    smooth = cv2.GaussianBlur(work, (0, 0), 12)
    bg = float(np.percentile(smooth, 5))
    hi = float(np.percentile(smooth, 95))
    threshold = bg + 0.12 * (hi - bg)
    mask = (smooth > threshold).astype(np.uint8) * 255
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, np.ones((31, 31), np.uint8))
    count, _labels, stats, _centroids = cv2.connectedComponentsWithStats(mask)

    candidates: list[tuple[int, int, int, int, int]] = []
    image_area = height * width
    for idx in range(1, count):
        area = int(stats[idx, cv2.CC_STAT_AREA])
        x = int(stats[idx, cv2.CC_STAT_LEFT])
        y = int(stats[idx, cv2.CC_STAT_TOP])
        w = int(stats[idx, cv2.CC_STAT_WIDTH])
        h = int(stats[idx, cv2.CC_STAT_HEIGHT])
        if area < image_area * 0.02:
            continue
        if h < height * 0.45 or y + h < height * 0.75:
            continue
        aspect = w / max(h, 1)
        if not 0.65 <= aspect <= 2.8:
            continue
        candidates.append((area, x, y, w, h))

    if not candidates:
        return PlateBox(0, int(height * 0.25), int(width * 0.80), int(height * 0.70))

    selected = sorted(candidates, key=lambda item: item[0], reverse=True)[:2]
    pad = 20
    left = max(0, min(x for _area, x, _y, _w, _h in selected) - pad)
    top = max(0, min(y for _area, _x, y, _w, _h in selected) - pad)
    right = min(width, max(x + w for _area, x, _y, w, _h in selected) + pad)
    bottom = min(height, max(y + h for _area, _x, y, _w, h in selected) + pad)
    return PlateBox(left, top, right - left, bottom - top)


def detect_dual_layer_centers(
    image: np.ndarray,
    roi_size: int,
    reverse_columns: bool = False,
) -> tuple[dict[str, dict[int, tuple[float, float]]], dict[str, PlateBox]]:
    height, width = image.shape
    array_box = detect_dual_layer_array_box(image)
    bottom_start = array_box.top + int(array_box.height * 0.55)
    bottom_stop = array_box.top + int(array_box.height * 0.98)
    bottom_band = image[
        bottom_start:bottom_stop,
        array_box.left:array_box.right,
    ].astype(np.float32)
    x_peaks = find_profile_peaks(
        bottom_band.mean(axis=0),
        expected=8,
        min_distance=max(int(array_box.width / 18), roi_size * 2),
        prominence_scale=0.12,
    )
    if len(x_peaks) != 8:
        raise RuntimeError(
            "Could not detect the two sMBC layer column grids; "
            "try manual processing or provide centers with --centers-json"
        )
    x_peaks = sorted(array_box.left + peak for peak in x_peaks)
    layer_xs = {
        "layer1": [float(x) for x in x_peaks[:4]],
        "layer2": [float(x) for x in x_peaks[4:8]],
    }

    y_band = image[
        array_box.top:array_box.bottom,
        array_box.left:array_box.right,
    ].astype(np.float32)
    y_peaks = find_profile_peaks(
        y_band.mean(axis=1),
        expected=4,
        min_distance=max(int(array_box.height / 7), roi_size * 2),
        prominence_scale=0.12,
    )
    if len(y_peaks) != 4:
        ys_bottom_up = fallback_grid_from_box(array_box)[1]
    else:
        ys_bottom_up = sorted((float(array_box.top + y) for y in y_peaks), reverse=True)

    centers_by_layer: dict[str, dict[int, tuple[float, float]]] = {}
    boxes_by_layer: dict[str, PlateBox] = {}
    for layer, xs in layer_xs.items():
        row_xs = list(reversed(xs)) if reverse_columns else xs
        channel_centers: dict[int, tuple[float, float]] = {}
        channel = 1
        for y in ys_bottom_up:
            for x in row_xs:
                channel_centers[channel] = refine_center(image, (x, y), roi_size)
                channel += 1
        centers_by_layer[layer] = channel_centers
        boxes_by_layer[layer] = box_from_grid(xs, ys_bottom_up, image.shape)

    return centers_by_layer, boxes_by_layer


def detect_layer_centers(
    image: np.ndarray,
    shot: str,
    roi_size: int,
    reverse_columns: bool = False,
) -> tuple[dict[str, dict[int, tuple[float, float]]], dict[str, PlateBox]]:
    if is_dual_layer_shot(shot):
        return detect_dual_layer_centers(image, roi_size, reverse_columns)
    centers, box = detect_centers(image, roi_size, reverse_columns)
    return {"layer1": centers}, {"layer1": box}


def parse_channel_centers(raw: dict[str, object], path: Path) -> dict[int, tuple[float, float]]:
    centers: dict[int, tuple[float, float]] = {}
    for key, value in raw.items():
        channel = int(key)
        if channel not in CHANNELS:
            raise ValueError(f"Invalid channel in {path}: {key}")
        if not isinstance(value, dict):
            raise ValueError(f"Invalid center entry in {path}: {key}")
        centers[channel] = (float(value["x"]), float(value["y"]))
    if sorted(centers) != CHANNELS:
        raise ValueError(f"{path} must contain channels 1..16")
    return centers


def load_centers_json(path: Path) -> dict[str, dict[int, tuple[float, float]]]:
    raw = json.loads(path.read_text())
    if all(str(channel) in raw for channel in CHANNELS):
        return {"layer1": parse_channel_centers(raw, path)}
    centers_by_layer: dict[str, dict[int, tuple[float, float]]] = {}
    for layer, centers in raw.items():
        if not isinstance(centers, dict):
            raise ValueError(f"Invalid layer entry in {path}: {layer}")
        centers_by_layer[str(layer)] = parse_channel_centers(centers, path)
    if not centers_by_layer:
        raise ValueError(f"No centers found in {path}")
    return centers_by_layer


def save_centers_json(
    path: Path,
    centers_by_layer: dict[str, dict[int, tuple[float, float]]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if list(centers_by_layer) == ["layer1"]:
        centers = centers_by_layer["layer1"]
        raw = {
            str(channel): {"x": float(centers[channel][0]), "y": float(centers[channel][1])}
            for channel in CHANNELS
        }
    else:
        raw = {}
        for layer, centers in centers_by_layer.items():
            raw[layer] = {
                str(channel): {"x": float(centers[channel][0]), "y": float(centers[channel][1])}
                for channel in CHANNELS
            }
    path.write_text(json.dumps(raw, indent=2))


class CenterEditor:
    def __init__(
        self,
        image: np.ndarray,
        centers_by_layer: dict[str, dict[int, tuple[float, float]]],
        roi_size: int,
        shot: str,
        boxes_by_layer: dict[str, PlateBox],
        cmap: str,
    ) -> None:
        self.image = image
        self.centers_by_layer = {
            layer: dict(centers)
            for layer, centers in sorted(centers_by_layer.items(), key=lambda item: layer_sort_key(item[0]))
        }
        self.roi_size = roi_size
        self.shot = shot
        self.boxes_by_layer = boxes_by_layer
        self.view_box = union_plate_box(list(boxes_by_layer.values()))
        self.cmap = cmap
        self.selected: tuple[str, int] | None = None
        self.drag_start: tuple[float, float] | None = None
        self.original_centers: dict[str, dict[int, tuple[float, float]]] | None = None
        self.move_all = False
        self.accepted = False
        self.fig, self.ax = plt.subplots(figsize=(10, 10))
        self.artists: dict[tuple[str, int], tuple[Rectangle, object, object]] = {}

    def run(self) -> dict[str, dict[int, tuple[float, float]]]:
        lo, hi = display_limits(self.image)
        image_artist = self.ax.imshow(
            self.image,
            cmap=self.cmap,
            vmin=lo,
            vmax=hi,
            origin="upper",
        )
        colorbar = self.fig.colorbar(image_artist, ax=self.ax, fraction=0.046, pad=0.04)
        colorbar.set_label("intensity")
        self.ax.set_title(
            "sMBC shot "
            + self.shot
            + " | drag labels/boxes; m toggles move-all; Enter accepts; q abort"
        )
        set_plate_view(self.ax, self.image.shape, self.view_box)
        for layer, box in self.boxes_by_layer.items():
            self.ax.add_patch(
                Rectangle(
                    (box.left, box.top),
                    box.width,
                    box.height,
                    fill=False,
                    edgecolor="cyan" if layer == "layer1" else "magenta",
                    linewidth=1.0,
                    linestyle="--",
                )
            )
            self.ax.text(box.left, max(0, box.top - 12), layer, color="white", fontsize=10)
        self.draw_centers()
        self.fig.canvas.mpl_connect("button_press_event", self.on_press)
        self.fig.canvas.mpl_connect("button_release_event", self.on_release)
        self.fig.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        plt.show()
        if not self.accepted:
            raise RuntimeError("Manual adjustment aborted")
        return self.centers_by_layer

    def draw_centers(self) -> None:
        for artists in self.artists.values():
            for artist in artists:
                artist.remove()
        self.artists.clear()
        half = self.roi_size / 2
        for layer, centers in self.centers_by_layer.items():
            for channel in CHANNELS:
                key = (layer, channel)
                x, y = centers[channel]
                color = "red" if key == self.selected else "yellow"
                rect = Rectangle(
                    (x - half, y - half),
                    self.roi_size,
                    self.roi_size,
                    fill=False,
                    edgecolor=color,
                    linewidth=1.7,
                )
                point = self.ax.plot(x, y, marker="+", color=color, markersize=9)[0]
                label = self.ax.text(
                    x + half + 4,
                    y,
                    channel_label(layer, channel),
                    color=color,
                    fontsize=10,
                    weight="bold",
                    va="center",
                )
                self.ax.add_patch(rect)
                self.artists[key] = (rect, point, label)
        self.fig.canvas.draw_idle()

    def nearest_channel(self, x: float, y: float) -> tuple[str, int] | None:
        distances: list[tuple[float, tuple[str, int]]] = []
        for layer, centers in self.centers_by_layer.items():
            for channel in CHANNELS:
                distance = math.hypot(centers[channel][0] - x, centers[channel][1] - y)
                distances.append((distance, (layer, channel)))
        distance, key = min(distances)
        return key if distance <= max(self.roi_size * 2, 80) else None

    def on_press(self, event) -> None:  # noqa: ANN001
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return
        key = self.nearest_channel(float(event.xdata), float(event.ydata))
        if key is None:
            return
        self.selected = key
        self.drag_start = (float(event.xdata), float(event.ydata))
        self.original_centers = {
            layer: dict(centers)
            for layer, centers in self.centers_by_layer.items()
        }
        self.draw_centers()

    def on_release(self, _event) -> None:  # noqa: ANN001
        self.drag_start = None
        self.original_centers = None

    def on_motion(self, event) -> None:  # noqa: ANN001
        if (
            self.selected is None
            or self.drag_start is None
            or self.original_centers is None
            or event.inaxes != self.ax
            or event.xdata is None
            or event.ydata is None
        ):
            return
        dx = float(event.xdata) - self.drag_start[0]
        dy = float(event.ydata) - self.drag_start[1]
        if self.move_all:
            for layer, centers in self.original_centers.items():
                for channel in CHANNELS:
                    ox, oy = centers[channel]
                    self.centers_by_layer[layer][channel] = (ox + dx, oy + dy)
        else:
            layer, channel = self.selected
            ox, oy = self.original_centers[layer][channel]
            self.centers_by_layer[layer][channel] = (ox + dx, oy + dy)
        self.draw_centers()

    def on_key(self, event) -> None:  # noqa: ANN001
        if event.key in ("enter", "return"):
            self.accepted = True
            plt.close(self.fig)
            return
        if event.key in ("q", "escape"):
            self.accepted = False
            plt.close(self.fig)
            return
        if event.key == "m":
            self.move_all = not self.move_all
            self.ax.set_title(
                "sMBC shot "
                + self.shot
                + (" | MOVE-ALL enabled" if self.move_all else " | move one channel")
                + " | Enter accepts; q abort"
            )
            self.fig.canvas.draw_idle()
            return
        if self.selected is None:
            return
        step = 10 if str(event.key).startswith("shift+") else 1
        key = str(event.key).replace("shift+", "")
        moves = {
            "left": (-step, 0),
            "right": (step, 0),
            "up": (0, -step),
            "down": (0, step),
        }
        if key not in moves:
            return
        dx, dy = moves[key]
        if self.move_all:
            selected_items = [
                (layer, channel)
                for layer in self.centers_by_layer
                for channel in CHANNELS
            ]
        else:
            selected_items = [self.selected]
        for layer, channel in selected_items:
            x, y = self.centers_by_layer[layer][channel]
            self.centers_by_layer[layer][channel] = (x + dx, y + dy)
        self.draw_centers()


def extract_intensities(
    image: np.ndarray,
    centers: dict[int, tuple[float, float]],
    roi_size: int,
    layer: str = "layer1",
) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    half = roi_size // 2
    for channel in CHANNELS:
        x, y = centers[channel]
        cx = int(round(x))
        cy = int(round(y))
        x1 = max(0, cx - half)
        y1 = max(0, cy - half)
        x2 = min(image.shape[1], x1 + roi_size)
        y2 = min(image.shape[0], y1 + roi_size)
        if x2 - x1 < roi_size:
            x1 = max(0, x2 - roi_size)
        if y2 - y1 < roi_size:
            y1 = max(0, y2 - roi_size)
        patch = image[y1:y2, x1:x2].astype(np.float64)
        rows.append(
            {
                "layer": layer,
                "channel": channel,
                "center_x": float(x),
                "center_y": float(y),
                "roi_x1": int(x1),
                "roi_y1": int(y1),
                "roi_x2": int(x2),
                "roi_y2": int(y2),
                "intensity_mean": float(np.mean(patch)),
                "intensity_std": float(np.std(patch, ddof=1)) if patch.size > 1 else 0.0,
                "intensity_min": float(np.min(patch)),
                "intensity_max": float(np.max(patch)),
            }
        )
    return pd.DataFrame(rows)


def extract_layer_intensities(
    image: np.ndarray,
    centers_by_layer: dict[str, dict[int, tuple[float, float]]],
    roi_size: int,
) -> pd.DataFrame:
    frames = [
        extract_intensities(image, centers, roi_size, layer)
        for layer, centers in sorted(centers_by_layer.items(), key=lambda item: layer_sort_key(item[0]))
    ]
    return pd.concat(frames, ignore_index=True)


def write_excel(
    excel_path: Path,
    shot: str,
    measurements: pd.DataFrame,
    image_info: ImageInfo,
    roi_size: int,
) -> None:
    excel_path.parent.mkdir(parents=True, exist_ok=True)
    now = dt.datetime.now().isoformat(timespec="seconds")
    long_df = measurements.copy()
    long_df.insert(0, "shot", str(shot))
    if "layer" not in long_df.columns:
        long_df.insert(1, "layer", "layer1")
    long_df["image_path"] = str(image_info.path)
    long_df["roi_size_px"] = roi_size
    long_df["img_dtype"] = image_info.dtype
    long_df["processed_at"] = now

    if excel_path.exists():
        try:
            existing_long = pd.read_excel(excel_path, sheet_name="intensity_long", dtype={"shot": str})
        except ValueError:
            existing_long = pd.DataFrame()
    else:
        existing_long = pd.DataFrame()

    if not existing_long.empty:
        if "layer" not in existing_long.columns:
            existing_long.insert(1, "layer", "layer1")
        existing_long = existing_long[existing_long["shot"].astype(str) != str(shot)]
    all_long = pd.concat([existing_long, long_df], ignore_index=True)
    all_long["channel"] = all_long["channel"].astype(int)
    all_long["layer"] = all_long["layer"].fillna("layer1").astype(str)
    all_long["_layer_order"] = all_long["layer"].map(lambda value: layer_sort_key(str(value))[0])
    all_long = all_long.sort_values(["shot", "_layer_order", "channel"]).drop(columns="_layer_order").reset_index(drop=True)

    wide_rows: list[dict[str, object]] = []
    for shot_value, group in all_long.groupby("shot", sort=True):
        row: dict[str, object] = {"shot": shot_value}
        for _idx, measurement in group.iterrows():
            layer = str(measurement["layer"])
            channel = int(measurement["channel"])
            column = f"ch{channel:02d}" if layer == "layer1" else f"ch{channel:02d}_{layer.removeprefix('layer')}"
            row[column] = float(measurement["intensity_mean"])
        wide_rows.append(row)
    wide = pd.DataFrame(wide_rows)
    ordered_columns = ["shot"]
    ordered_columns.extend(f"ch{channel:02d}" for channel in CHANNELS)
    ordered_columns.extend(f"ch{channel:02d}_2" for channel in CHANNELS)
    ordered_columns.extend(column for column in wide.columns if column not in ordered_columns)
    wide = wide.reindex(columns=[column for column in ordered_columns if column in wide.columns])

    processed_layers = ",".join(
        sorted(long_df["layer"].astype(str).unique(), key=lambda value: layer_sort_key(value)[0])
    )

    metadata = pd.DataFrame(
        [
            {
                "shot": str(shot),
                "image_path": str(image_info.path),
                "inf_path": str(image_info.inf_path) if image_info.inf_path else "",
                "width": image_info.width,
                "height": image_info.height,
                "img_dtype": image_info.dtype,
                "roi_size_px": roi_size,
                "layers": processed_layers,
                "processed_at": now,
            }
        ]
    )
    if excel_path.exists():
        try:
            existing_meta = pd.read_excel(excel_path, sheet_name="metadata", dtype={"shot": str})
            existing_meta = existing_meta[existing_meta["shot"].astype(str) != str(shot)]
            metadata = pd.concat([existing_meta, metadata], ignore_index=True)
        except ValueError:
            pass
    metadata = metadata.sort_values("shot").reset_index(drop=True)

    with pd.ExcelWriter(excel_path, engine="openpyxl", mode="w") as writer:
        all_long.to_excel(writer, index=False, sheet_name="intensity_long")
        wide.to_excel(writer, index=False, sheet_name="intensity_wide")
        metadata.to_excel(writer, index=False, sheet_name="metadata")


def save_preview(
    output_path: Path,
    image: np.ndarray,
    centers_by_layer: dict[str, dict[int, tuple[float, float]]],
    roi_size: int,
    shot: str,
    boxes_by_layer: dict[str, PlateBox] | None = None,
    cmap: str = "jet",
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lo, hi = display_limits(image)
    fig, ax = plt.subplots(figsize=(10, 10))
    image_artist = ax.imshow(image, cmap=cmap, vmin=lo, vmax=hi, origin="upper")
    colorbar = fig.colorbar(image_artist, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("intensity")
    ax.set_title(f"sMBC shot {shot} center preview")
    ax.axis("off")
    if boxes_by_layer is not None:
        set_plate_view(ax, image.shape, union_plate_box(list(boxes_by_layer.values())))
        for layer, box in boxes_by_layer.items():
            ax.add_patch(
                Rectangle(
                    (box.left, box.top),
                    box.width,
                    box.height,
                    fill=False,
                    edgecolor="cyan" if layer == "layer1" else "magenta",
                    linewidth=1.0,
                    linestyle="--",
                )
            )
            ax.text(box.left, max(0, box.top - 12), layer, color="white", fontsize=10)
    half = roi_size / 2
    for layer, centers in sorted(centers_by_layer.items(), key=lambda item: layer_sort_key(item[0])):
        for channel in CHANNELS:
            x, y = centers[channel]
            ax.add_patch(
                Rectangle(
                    (x - half, y - half),
                    roi_size,
                    roi_size,
                    fill=False,
                    edgecolor="yellow",
                    linewidth=1.6,
                )
            )
            ax.plot(x, y, marker="+", color="red", markersize=8)
            ax.text(
                x + half + 4,
                y,
                channel_label(layer, channel),
                color="yellow",
                fontsize=10,
                weight="bold",
            )
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def process_shot(args: argparse.Namespace) -> None:
    image_path = locate_image(args.shot, args.data_root, args.image)
    image, image_info = load_img(image_path, args.byte_order)
    if args.centers_json:
        centers_by_layer = load_centers_json(args.centers_json)
        boxes_by_layer = boxes_from_centers(centers_by_layer, image.shape)
    else:
        centers_by_layer, boxes_by_layer = detect_layer_centers(
            image,
            args.shot,
            args.roi_size,
            args.reverse_columns,
        )

    if not args.no_gui:
        centers_by_layer = CenterEditor(
            image,
            centers_by_layer,
            args.roi_size,
            args.shot,
            boxes_by_layer,
            args.preview_cmap,
        ).run()
        boxes_by_layer = boxes_from_centers(centers_by_layer, image.shape)

    measurements = extract_layer_intensities(image, centers_by_layer, args.roi_size)
    write_excel(args.excel, args.shot, measurements, image_info, args.roi_size)

    preview = args.preview or args.excel.with_name(f"sMBC_{args.shot}_preview.png")
    save_preview(
        preview,
        image,
        centers_by_layer,
        args.roi_size,
        args.shot,
        boxes_by_layer,
        args.preview_cmap,
    )
    if args.save_centers_json:
        save_centers_json(args.save_centers_json, centers_by_layer)

    print(f"Processed shot {args.shot}")
    print(f"Image: {image_path}")
    print(f"Excel: {args.excel}")
    print(f"Preview: {preview}")
    print(measurements[["layer", "channel", "center_x", "center_y", "intensity_mean"]].to_string(index=False))


def normalize_long_dataframe(long_df: pd.DataFrame) -> pd.DataFrame:
    normalized = long_df.copy()
    if "layer" not in normalized.columns:
        normalized.insert(1, "layer", "layer1")
    normalized["shot"] = normalized["shot"].astype(str)
    normalized["layer"] = normalized["layer"].fillna("layer1").astype(str)
    normalized["channel"] = normalized["channel"].astype(int)
    return normalized


def resolve_plot_shots(args: argparse.Namespace, long_df: pd.DataFrame) -> tuple[list[str], list[str]]:
    if args.shot_range and args.shots:
        raise RuntimeError("Use either explicit shots or --range, not both")
    if not args.shot_range and not args.shots:
        raise RuntimeError("Provide at least one shot or use --range START END")

    available = set(long_df["shot"].astype(str))
    if args.shot_range:
        start, end = (int(value) for value in args.shot_range)
        if start > end:
            raise RuntimeError("--range START END requires START <= END")
        requested = [str(shot) for shot in range(start, end + 1)]
        missing = [shot for shot in requested if shot not in available]
        if missing:
            print(
                "Warning: skipping shots not found in Excel: " + ", ".join(missing),
                file=sys.stderr,
            )
        shots = [shot for shot in requested if shot in available]
        descriptor = [str(start), str(end)]
    else:
        shots = [str(shot) for shot in args.shots]
        missing = [shot for shot in shots if shot not in available]
        if missing:
            raise RuntimeError(f"Shots not found in {args.excel}: {', '.join(missing)}")
        descriptor = shots

    if not shots:
        raise RuntimeError("No matching shots found to plot")
    return shots, descriptor


def resolve_plot_groups(group_args: list[str], shots: list[str]) -> list[list[str]]:
    shot_set = set(shots)
    assigned: dict[str, int] = {}
    groups: list[list[str]] = []

    for group_index, raw_group in enumerate(group_args, start=1):
        group_shots = [shot.strip() for shot in raw_group.split(",") if shot.strip()]
        if not group_shots:
            raise RuntimeError(f"--group #{group_index} is empty")

        valid_group: list[str] = []
        missing: list[str] = []
        for shot in group_shots:
            if shot not in shot_set:
                missing.append(shot)
                continue
            if shot in assigned:
                raise RuntimeError(
                    f"Shot {shot} appears in multiple groups "
                    f"(groups {assigned[shot]} and {group_index})"
                )
            assigned[shot] = group_index
            valid_group.append(shot)

        if missing:
            print(
                f"Warning: ignoring --group #{group_index} shots not in plotted set: "
                + ", ".join(missing),
                file=sys.stderr,
            )
        if valid_group:
            groups.append(valid_group)

    for shot in shots:
        if shot not in assigned:
            groups.append([shot])

    return groups


def plot_shots(args: argparse.Namespace) -> None:
    if not args.excel.exists():
        raise FileNotFoundError(f"Excel file not found: {args.excel}")
    long_df = normalize_long_dataframe(
        pd.read_excel(args.excel, sheet_name="intensity_long", dtype={"shot": str})
    )
    shots, descriptor = resolve_plot_shots(args, long_df)
    subset = long_df[(long_df["shot"].isin(shots)) & (long_df["layer"] == args.layer)].copy()
    missing_layer = [
        shot
        for shot in shots
        if subset[subset["shot"].astype(str) == shot].empty
    ]
    if missing_layer:
        print(
            f"Warning: skipping shots without {args.layer} data: " + ", ".join(missing_layer),
            file=sys.stderr,
        )
        shots = [shot for shot in shots if shot not in missing_layer]
    if not shots:
        raise RuntimeError(f"No {args.layer} data found to plot")

    plot_groups = resolve_plot_groups(args.group, shots)
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    marker_cycle = ["o", "s", "^", "D", "v", "P", "X", "*"]

    fig, ax = plt.subplots(figsize=(8, 5))
    for group_index, group_shots in enumerate(plot_groups):
        color = color_cycle[group_index % len(color_cycle)] if color_cycle else None
        for shot_index, shot in enumerate(group_shots):
            data = subset[subset["shot"].astype(str) == shot].sort_values("channel")
            ax.plot(
                data["channel"],
                data["intensity_mean"],
                marker=marker_cycle[shot_index % len(marker_cycle)],
                color=color,
                linewidth=1.8,
                label=shot,
            )
    ax.set_xlabel("channel")
    ax.set_ylabel(f"{args.layer} intensity mean (50x50 px)")
    if args.log_y:
        if (subset["intensity_mean"] <= 0).any():
            raise RuntimeError("Cannot use --log-y because at least one intensity is <= 0")
        ax.set_yscale("log")
    ax.set_xticks(CHANNELS)
    ax.grid(True, which="both" if args.log_y else "major", alpha=0.25)
    ax.legend(title=f"{args.layer} shot")
    fig.tight_layout()

    output = args.output or args.excel.with_name(
        compare_plot_name(descriptor, args.layer, args.log_y, bool(args.group))
    )
    output = make_unique_path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    if args.show:
        plt.show()
    else:
        plt.close(fig)
    print(f"Plot saved: {output}")


def shot_sort_key(shot: str) -> tuple[int, int, str]:
    shot_number = parse_shot_number(shot)
    if shot_number is None:
        return (1, 0, str(shot))
    return (0, shot_number, str(shot))


def read_available_plot_shots(excel_path: Path) -> list[str]:
    if not excel_path.exists():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    long_df = normalize_long_dataframe(
        pd.read_excel(excel_path, sheet_name="intensity_long", dtype={"shot": str})
    )
    return sorted(long_df["shot"].astype(str).unique().tolist(), key=shot_sort_key)


class PlotGui:
    def __init__(self, initial_excel: Path) -> None:
        import tkinter as tk
        from tkinter import filedialog, scrolledtext, ttk

        self.tk = tk
        self.filedialog = filedialog
        self.scrolledtext = scrolledtext
        self.ttk = ttk
        self.available_shots: list[str] = []
        self.groups: list[list[str]] = []

        self.root = tk.Tk()
        self.root.title("sMBC Plot GUI")
        self.root.minsize(780, 560)

        self.excel_var = tk.StringVar(value=str(initial_excel))
        self.output_var = tk.StringVar(value="")
        self.layer_var = tk.StringVar(value="layer1")
        self.log_y_var = tk.BooleanVar(value=False)
        self.show_var = tk.BooleanVar(value=True)

        self.build_layout()
        self.load_excel()

    def run(self) -> None:
        self.root.mainloop()

    def build_layout(self) -> None:
        tk = self.tk
        ttk = self.ttk

        outer = ttk.Frame(self.root, padding=12)
        outer.grid(row=0, column=0, sticky="nsew")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(1, weight=1)
        outer.rowconfigure(4, weight=1)

        excel_frame = ttk.Frame(outer)
        excel_frame.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        excel_frame.columnconfigure(1, weight=1)
        ttk.Label(excel_frame, text="Excel").grid(row=0, column=0, sticky="w")
        ttk.Entry(excel_frame, textvariable=self.excel_var).grid(
            row=0,
            column=1,
            sticky="ew",
            padx=8,
        )
        ttk.Button(excel_frame, text="Browse", command=self.browse_excel).grid(
            row=0,
            column=2,
            padx=(0, 6),
        )
        ttk.Button(excel_frame, text="Load", command=self.load_excel).grid(row=0, column=3)

        selector = ttk.Frame(outer)
        selector.grid(row=1, column=0, sticky="nsew")
        selector.columnconfigure(0, weight=1)
        selector.columnconfigure(2, weight=1)
        selector.rowconfigure(1, weight=1)

        ttk.Label(selector, text="Available shots").grid(row=0, column=0, sticky="w")
        ttk.Label(selector, text="Groups").grid(row=0, column=2, sticky="w")

        shot_frame = ttk.Frame(selector)
        shot_frame.grid(row=1, column=0, sticky="nsew")
        shot_frame.columnconfigure(0, weight=1)
        shot_frame.rowconfigure(0, weight=1)
        self.available_list = tk.Listbox(
            shot_frame,
            selectmode=tk.EXTENDED,
            exportselection=False,
            height=14,
        )
        self.available_list.grid(row=0, column=0, sticky="nsew")
        self.available_list.bind("<Double-Button-1>", lambda _event: self.add_group())
        shot_scroll = ttk.Scrollbar(
            shot_frame,
            orient=tk.VERTICAL,
            command=self.available_list.yview,
        )
        shot_scroll.grid(row=0, column=1, sticky="ns")
        self.available_list.configure(yscrollcommand=shot_scroll.set)

        add_frame = ttk.Frame(selector)
        add_frame.grid(row=1, column=1, sticky="n", padx=12)
        ttk.Button(add_frame, text="Add group >", command=self.add_group).grid(
            row=0,
            column=0,
            sticky="ew",
            pady=(28, 0),
        )

        group_frame = ttk.Frame(selector)
        group_frame.grid(row=1, column=2, sticky="nsew")
        group_frame.columnconfigure(0, weight=1)
        group_frame.rowconfigure(0, weight=1)
        self.group_list = tk.Listbox(
            group_frame,
            selectmode=tk.EXTENDED,
            exportselection=False,
            height=14,
        )
        self.group_list.grid(row=0, column=0, sticky="nsew")
        group_scroll = ttk.Scrollbar(
            group_frame,
            orient=tk.VERTICAL,
            command=self.group_list.yview,
        )
        group_scroll.grid(row=0, column=1, sticky="ns")
        self.group_list.configure(yscrollcommand=group_scroll.set)

        group_buttons = ttk.Frame(group_frame)
        group_buttons.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        for column in range(4):
            group_buttons.columnconfigure(column, weight=1)
        ttk.Button(group_buttons, text="Remove", command=self.remove_groups).grid(
            row=0,
            column=0,
            sticky="ew",
            padx=(0, 4),
        )
        ttk.Button(group_buttons, text="Clear", command=self.clear_groups).grid(
            row=0,
            column=1,
            sticky="ew",
            padx=4,
        )
        ttk.Button(group_buttons, text="Up", command=lambda: self.move_group(-1)).grid(
            row=0,
            column=2,
            sticky="ew",
            padx=4,
        )
        ttk.Button(group_buttons, text="Down", command=lambda: self.move_group(1)).grid(
            row=0,
            column=3,
            sticky="ew",
            padx=(4, 0),
        )

        options = ttk.Frame(outer)
        options.grid(row=2, column=0, sticky="ew", pady=(12, 8))
        options.columnconfigure(5, weight=1)
        ttk.Label(options, text="Layer").grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(
            options,
            text="layer1",
            value="layer1",
            variable=self.layer_var,
        ).grid(row=0, column=1, sticky="w", padx=(8, 0))
        ttk.Radiobutton(
            options,
            text="layer2",
            value="layer2",
            variable=self.layer_var,
        ).grid(row=0, column=2, sticky="w", padx=(8, 0))
        ttk.Checkbutton(options, text="log-y", variable=self.log_y_var).grid(
            row=0,
            column=3,
            sticky="w",
            padx=(16, 0),
        )
        ttk.Checkbutton(options, text="show plot", variable=self.show_var).grid(
            row=0,
            column=4,
            sticky="w",
            padx=(16, 0),
        )

        output_frame = ttk.Frame(outer)
        output_frame.grid(row=3, column=0, sticky="ew", pady=(0, 10))
        output_frame.columnconfigure(1, weight=1)
        ttk.Label(output_frame, text="Output").grid(row=0, column=0, sticky="w")
        ttk.Entry(output_frame, textvariable=self.output_var).grid(
            row=0,
            column=1,
            sticky="ew",
            padx=8,
        )
        ttk.Button(output_frame, text="Browse", command=self.browse_output).grid(
            row=0,
            column=2,
        )

        status_frame = ttk.Frame(outer)
        status_frame.grid(row=4, column=0, sticky="nsew")
        status_frame.columnconfigure(0, weight=1)
        status_frame.rowconfigure(1, weight=1)
        ttk.Button(status_frame, text="Run plot", command=self.run_plot).grid(
            row=0,
            column=0,
            sticky="e",
            pady=(0, 8),
        )
        self.status_text = self.scrolledtext.ScrolledText(
            status_frame,
            height=8,
            wrap=tk.WORD,
            state="disabled",
        )
        self.status_text.grid(row=1, column=0, sticky="nsew")

    def browse_excel(self) -> None:
        initial_path = Path(self.excel_var.get()).expanduser()
        initial_dir = initial_path.parent if initial_path.parent.exists() else Path.cwd()
        selected = self.filedialog.askopenfilename(
            title="Select sMBC intensity Excel file",
            initialdir=str(initial_dir),
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")],
        )
        if selected:
            self.excel_var.set(selected)
            self.load_excel()

    def browse_output(self) -> None:
        excel_path = Path(self.excel_var.get()).expanduser()
        initial_dir = excel_path.parent if excel_path.parent.exists() else Path.cwd()
        selected = self.filedialog.asksaveasfilename(
            title="Save comparison plot",
            initialdir=str(initial_dir),
            initialfile="sMBC_compare.png",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")],
        )
        if selected:
            self.output_var.set(selected)

    def load_excel(self) -> None:
        excel_path = Path(self.excel_var.get()).expanduser()
        try:
            shots = read_available_plot_shots(excel_path)
        except Exception as exc:  # noqa: BLE001
            self.available_shots = []
            self.groups = []
            self.refresh_available_shots()
            self.refresh_groups()
            self.set_status(f"Error loading Excel: {exc}")
            return

        self.available_shots = shots
        self.groups = []
        self.refresh_available_shots()
        self.refresh_groups()
        self.set_status(f"Loaded {len(shots)} shots from {excel_path}")

    def refresh_available_shots(self) -> None:
        self.available_list.delete(0, self.tk.END)
        for shot in self.available_shots:
            self.available_list.insert(self.tk.END, shot)

    def refresh_groups(self) -> None:
        self.group_list.delete(0, self.tk.END)
        for index, group in enumerate(self.groups, start=1):
            self.group_list.insert(self.tk.END, f"Group {index}: {', '.join(group)}")

    def grouped_shots(self) -> set[str]:
        return {shot for group in self.groups for shot in group}

    def add_group(self) -> None:
        selected_indices = self.available_list.curselection()
        if not selected_indices:
            self.set_status("Select one or more shots before adding a group.")
            return

        selected = [self.available_shots[index] for index in selected_indices]
        already_grouped = [shot for shot in selected if shot in self.grouped_shots()]
        if already_grouped:
            self.set_status(
                "These shots are already in a group: " + ", ".join(already_grouped)
            )
            return

        self.groups.append(selected)
        self.refresh_groups()
        self.set_status("Added group: " + ", ".join(selected))

    def remove_groups(self) -> None:
        selected_indices = list(self.group_list.curselection())
        if not selected_indices:
            self.set_status("Select one or more groups to remove.")
            return
        for index in sorted(selected_indices, reverse=True):
            del self.groups[index]
        self.refresh_groups()
        self.set_status("Removed selected groups.")

    def clear_groups(self) -> None:
        self.groups = []
        self.refresh_groups()
        self.set_status("Cleared all groups.")

    def move_group(self, direction: int) -> None:
        selected_indices = list(self.group_list.curselection())
        if len(selected_indices) != 1:
            self.set_status("Select exactly one group to move.")
            return
        index = selected_indices[0]
        new_index = index + direction
        if new_index < 0 or new_index >= len(self.groups):
            return
        self.groups[index], self.groups[new_index] = self.groups[new_index], self.groups[index]
        self.refresh_groups()
        self.group_list.selection_set(new_index)
        self.group_list.activate(new_index)
        self.set_status("Moved selected group.")

    def build_plot_args(self) -> argparse.Namespace:
        if not self.groups:
            raise RuntimeError("Add at least one group before plotting")

        shots: list[str] = []
        group_args: list[str] = []
        for group in self.groups:
            if not group:
                continue
            group_args.append(",".join(group))
            for shot in group:
                if shot not in shots:
                    shots.append(shot)
        if not shots:
            raise RuntimeError("Add at least one shot before plotting")

        output_text = self.output_var.get().strip()
        output_path = Path(output_text).expanduser() if output_text else None

        return argparse.Namespace(
            command="plot",
            shots=shots,
            shot_range=None,
            excel=Path(self.excel_var.get()).expanduser(),
            output=output_path,
            show=bool(self.show_var.get()),
            log_y=bool(self.log_y_var.get()),
            group=group_args,
            layer=self.layer_var.get(),
        )

    def run_plot(self) -> None:
        try:
            plot_args = self.build_plot_args()
        except Exception as exc:  # noqa: BLE001
            self.set_status(f"Error: {exc}")
            return

        output = io.StringIO()
        try:
            self.set_status("Running plot...")
            self.root.update_idletasks()
            with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
                plot_shots(plot_args)
        except Exception as exc:  # noqa: BLE001
            details = output.getvalue().strip()
            message = f"Error: {exc}"
            if details:
                message += "\n\n" + details
            self.set_status(message)
            return

        message = output.getvalue().strip() or "Plot finished."
        self.set_status(message)

    def set_status(self, message: str) -> None:
        self.status_text.configure(state="normal")
        self.status_text.delete("1.0", self.tk.END)
        self.status_text.insert(self.tk.END, message.rstrip() + "\n")
        self.status_text.configure(state="disabled")


def run_plot_gui(args: argparse.Namespace) -> None:
    try:
        gui = PlotGui(args.excel)
    except ImportError as exc:
        raise RuntimeError("tkinter is not available in this Python environment") from exc
    gui.run()


def main() -> int:
    args = parse_args()
    try:
        if args.command == "process":
            process_shot(args)
        elif args.command == "plot":
            plot_shots(args)
        elif args.command == "gui-plot":
            run_plot_gui(args)
        else:
            raise AssertionError(args.command)
    except Exception as exc:  # noqa: BLE001
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
