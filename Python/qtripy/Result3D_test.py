"""Launch Result3D with a ventricular mesh and reproducible random values."""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
from dotenv import load_dotenv
from qtpy.QtWidgets import QApplication

PYTHON_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = PYTHON_ROOT.parent
sys.path.insert(0, str(PYTHON_ROOT))

from helpers.helpers import load_multi_array
from Result3D import Result3D


def _default_mesh_path() -> Path:
    load_dotenv(PYTHON_ROOT / ".env")
    data_root = Path(os.getenv("ENV_DATA_PATH", PROJECT_ROOT / "Data"))
    candidates = [
        data_root / "ECGsim_data" / "PPD1_ECGSIM_NEW" / "model" / "ventricle.tri",
        PROJECT_ROOT
        / "Data"
        / "raw"
        / "Models"
        / "ECGSIM_ALVALE_00_model"
        / "ALVALE_00_ventricles.tri",
    ]
    for path in candidates:
        if path.is_file():
            return path
    expected = "\n".join(f"  - {path}" for path in candidates)
    raise FileNotFoundError(
        "Could not find a ventricular mesh. Pass one with --mesh.\n"
        f"Checked:\n{expected}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mesh",
        type=Path,
        help="Path to a text .tri file containing vertex and triangle arrays",
    )
    parser.add_argument("--seed", type=int, default=2601)
    args = parser.parse_args()

    mesh_path = args.mesh or _default_mesh_path()
    arrays = load_multi_array(mesh_path)
    if len(arrays) < 2 or arrays[0].shape[1] < 4 or arrays[1].shape[1] < 4:
        raise ValueError(f"{mesh_path} must contain vertex and triangle blocks")

    vertices = arrays[0][:, 1:4]
    triangles = arrays[1][:, 1:4]
    random_values = np.arange(0, 100, 100/len(vertices))
    # random_values = np.random.default_rng(args.seed).uniform(0, 100, len(vertices))

    viewer = Result3D()
    viewer.begin()
    viewer.reset()
    viewer.surface(vertices, triangles, name="ventricle")
    viewer.transparency(0.8)
    viewer.values(random_values)
    viewer.gradient_bins()
    viewer.text("Ventricle - reproducible random values", (0.03, 0.95))
    viewer.background_color("white")
    viewer.property_on_mouse_click("coor")

    print(f"Displaying {len(vertices)} vertices from: {mesh_path}")
    print("Close the Result3D window to exit.")
    try:
        while viewer.plotter is not None and viewer.plotter.app_window.isVisible():
            QApplication.processEvents()
            time.sleep(0.05)
    finally:
        viewer.close()


if __name__ == "__main__":
    main()
