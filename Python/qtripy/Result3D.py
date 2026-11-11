"""Result3D: a PyVista-backed 3D results viewer.

The plotting methods are designed for MATLAB/Python workflows that display
triangulated geometry and scalar data in synchronized panels.

Install the viewer dependencies with ``pip install -r requirements-pyvista.txt``.

MATLAB can import this module after adding its directory to ``py.sys.path``::

    py.importlib.import_module('Result3D')
    viewer = py.Result3D.Result3D();
    viewer.begin();
"""

from __future__ import annotations

import sys
import time
from typing import Any

import numpy as np

if __package__:
    from .qtripy import QTripy
else:
    from qtripy import QTripy


class Result3D(QTripy):
    """Display triangulated models and scalar results in a PyVistaQt window."""

    def __init__(self, bins: int = 15) -> None:
        """Create a viewer with a default of 15 discrete scalar color bins."""
        self._validate_bins(bins)
        super().__init__(path=sys.executable)
        self.plotter: Any = None
        self._mesh_records: list[dict[str, Any]] = []
        self._text_records: list[dict[str, Any]] = []
        self._panel_count = (1, 1)
        self._background = "black"
        self._text_color = "black"
        self._edge_color: str | None = None
        self._opacity = 1.0
        self._color_map = "jet"
        self._color_limits: tuple[float, float] | None = None
        self._data_range: tuple[float, float] | None = None
        self._manual_color_limits = False
        self._bins = int(bins)
        self._lookup_table: Any = None
        self._scalar_bar_added = False
        self._controls_dock: Any = None
        self._range_min_field: Any = None
        self._range_max_field: Any = None
        self._bins_field: Any = None
        self._opacity_field: Any = None
        self._transparency_mode = "depth_peeling"
        self._controls_status: Any = None
        self._picked_node_number: Any = None
        self._picked_node_coordinates: Any = None
        self._picked_node_value: Any = None
        self._picked_node: tuple[dict[str, Any], int] | None = None
        self._updating_controls = False
        self._pick_property = "coor"
        self._frame_index = 0

    def begin(self, port: int = 1041) -> None:
        """Open the PyVista viewer; ``port`` is accepted for API compatibility."""
        del port
        if self.plotter is not None:
            return
        self._create_plotter()

    def _create_plotter(self) -> None:
        try:
            from pyvistaqt import BackgroundPlotter
        except ImportError as exc:
            raise ImportError(
                "Result3D requires PyVistaQt. Install dependencies with "
                "'pip install -r requirements-pyvista.txt'."
            ) from exc

        horizontal, vertical = self._panel_count
        self._scalar_bar_added = False
        self.plotter = BackgroundPlotter(
            title="Result3D",
            shape=(vertical, horizontal),
            auto_update=1,
        )
        self._create_controls_dock()
        self.plotter.set_background(self._background)
        self.plotter.add_axes()
        self.plotter.link_views()
        self._apply_transparency_mode()
        for record in self._mesh_records:
            self._draw_mesh(record)
        for record in self._text_records:
            self._draw_text(record)
        self._create_scalar_bar()
        self._update_controls()
        self.set_active_panel(*self.active_panel)

    @staticmethod
    def _make_polydata(vertices: Any, triangles: Any):
        import pyvista as pv

        points = np.asarray(vertices, dtype=float)
        faces_in = np.asarray(triangles, dtype=np.int64)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("vertices must have shape (N, 3)")
        if not np.isfinite(points).all():
            raise ValueError("vertices must contain only finite coordinates")
        if faces_in.size == 0:
            return pv.PolyData(points)
        if faces_in.ndim != 2 or faces_in.shape[1] != 3:
            raise ValueError("triangles must have shape (M, 3)")

        # MATLAB meshes are normally one-based, while Python-generated meshes
        # (such as QTripy's spheres) are zero-based.
        if faces_in.min() >= 1 and faces_in.max() <= len(points):
            faces_in = faces_in - 1
        if faces_in.min() < 0 or faces_in.max() >= len(points):
            raise ValueError("triangle indices are outside the vertices array")
        faces = np.column_stack(
            (np.full(len(faces_in), 3, dtype=np.int64), faces_in)
        ).ravel()
        return pv.PolyData(points, faces)

    def _ensure_started(self) -> None:
        if self.plotter is None:
            raise RuntimeError("Call begin() before sending plotting commands")

    def _should_debounce(self) -> bool:
        if self.min_time_before_next_cmd <= 0:
            return False
        current_time = time.time()
        if current_time - self.last_cmd_time < self.min_time_before_next_cmd:
            return True
        self.last_cmd_time = current_time
        return False

    def _select_panel(self, panel: tuple[int, int]) -> None:
        if self.plotter is not None:
            horizontal, vertical = panel
            self.plotter.subplot(vertical - 1, horizontal - 1)

    def _refresh_plotter(self) -> None:
        if self.plotter is None:
            return
        self._dispose_controls_dock()
        self.plotter.close()
        self.plotter = None
        self._create_plotter()

    def _create_controls_dock(self) -> None:
        from qtpy import QtCore, QtGui, QtWidgets

        self._controls_dock = QtWidgets.QDockWidget(
            "Result controls", self.plotter.app_window
        )
        self._controls_dock.setObjectName("Result3DControls")
        controls = QtWidgets.QWidget(self._controls_dock)
        layout = QtWidgets.QFormLayout(controls)

        self._range_min_field = QtWidgets.QLineEdit(controls)
        self._range_max_field = QtWidgets.QLineEdit(controls)
        self._bins_field = QtWidgets.QLineEdit(controls)
        self._opacity_field = QtWidgets.QLineEdit(controls)
        self._opacity_field.setToolTip(
            "1 is fully opaque and 0 is invisible. Values below 1 use "
            "depth-peeling transparency."
        )
        self._range_min_field.setValidator(
            QtGui.QDoubleValidator(-1e100, 1e100, 12, self._range_min_field)
        )
        self._range_max_field.setValidator(
            QtGui.QDoubleValidator(-1e100, 1e100, 12, self._range_max_field)
        )
        self._bins_field.setValidator(QtGui.QIntValidator(1, 1_000_000, self._bins_field))
        self._opacity_field.setValidator(
            QtGui.QDoubleValidator(0.0, 1.0, 6, self._opacity_field)
        )
        layout.addRow("Color minimum:", self._range_min_field)
        layout.addRow("Color maximum:", self._range_max_field)
        layout.addRow("Gradient bins:", self._bins_field)
        layout.addRow("Opacity (0-1):", self._opacity_field)

        apply_button = QtWidgets.QPushButton("Apply", controls)
        apply_button.clicked.connect(self._apply_display_settings)
        self._range_min_field.returnPressed.connect(self._apply_display_settings)
        self._range_max_field.returnPressed.connect(self._apply_display_settings)
        self._bins_field.returnPressed.connect(self._apply_display_settings)
        self._opacity_field.returnPressed.connect(self._apply_display_settings)
        layout.addRow(apply_button)

        self._controls_status = QtWidgets.QLabel(controls)
        self._controls_status.setWordWrap(True)
        layout.addRow(self._controls_status)

        node_group = QtWidgets.QGroupBox("Picked node", controls)
        node_layout = QtWidgets.QVBoxLayout(node_group)
        node_hint = QtWidgets.QLabel(
            "Right-click a model node to inspect it.",
            node_group,
        )
        node_layout.addWidget(node_hint)
        node_form = QtWidgets.QFormLayout()
        self._picked_node_number = QtWidgets.QLabel("—", node_group)
        self._picked_node_coordinates = QtWidgets.QLabel("—", node_group)
        self._picked_node_value = QtWidgets.QLabel("—", node_group)
        for label in (
            self._picked_node_number,
            self._picked_node_coordinates,
            self._picked_node_value,
        ):
            label.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
            label.setWordWrap(True)
        node_form.addRow("Node number:", self._picked_node_number)
        node_form.addRow("Coordinates:", self._picked_node_coordinates)
        node_form.addRow("Value:", self._picked_node_value)
        node_layout.addLayout(node_form)
        layout.addRow(node_group)
        self._reset_picked_node_info()

        self._controls_dock.setWidget(controls)
        self.plotter.app_window.addDockWidget(
            QtCore.Qt.DockWidgetArea.RightDockWidgetArea,
            self._controls_dock,
        )
        self._sync_display_fields()

    def _dispose_controls_dock(self) -> None:
        if self._controls_dock is None or self.plotter is None:
            return
        self.plotter.app_window.removeDockWidget(self._controls_dock)
        self._controls_dock.deleteLater()
        self._controls_dock = None
        self._range_min_field = None
        self._range_max_field = None
        self._bins_field = None
        self._opacity_field = None
        self._controls_status = None
        self._picked_node_number = None
        self._picked_node_coordinates = None
        self._picked_node_value = None

    def _sync_display_fields(self) -> None:
        if (
            self._range_min_field is None
            or self._range_max_field is None
            or self._opacity_field is None
        ):
            return
        limits = self._color_limits or self._data_range
        if limits is not None:
            self._range_min_field.setText(f"{limits[0]:.12g}")
            self._range_max_field.setText(f"{limits[1]:.12g}")
        if self._bins_field is not None:
            self._bins_field.setText(str(self._bins))
        self._opacity_field.setText(f"{self._opacity:.6g}")

    def _draw_mesh(self, record: dict[str, Any]) -> None:
        self._select_panel(record["panel"])
        mesh = record["mesh"]
        scalar_name = "values" if "values" in mesh.point_data else None
        options: dict[str, Any] = {
            "name": record["name"],
            "opacity": record["opacity"],
            "show_edges": record["show_edges"],
            "edge_color": self._normalize_color(self._edge_color or "black"),
            "cmap": self._color_map,
            "n_colors": self._bins,
            "show_scalar_bar": False,
        }
        if scalar_name is not None:
            options["scalars"] = scalar_name
            limits = record.get("color_limits") or self._color_limits
            if limits is not None:
                options["clim"] = limits
        else:
            options["color"] = record["color"]
        record["actor"] = self.plotter.add_mesh(mesh, **options)
        if scalar_name is not None and self._lookup_table is not None:
            record["actor"].mapper.lookup_table = self._lookup_table

    def _draw_text(self, record: dict[str, Any]) -> None:
        self._select_panel(record["panel"])
        record["actor"] = self.plotter.add_text(
            record["text"],
            position=record["position"],
            color=record["color"],
            font_size=record["font_size"],
        )

    def _find_mesh(self, name: str = "") -> dict[str, Any]:
        if not self._mesh_records:
            raise RuntimeError("No surface is available to receive scalar values")
        if name:
            exact = [record for record in self._mesh_records if record["name"] == name]
            matches = exact or [
                record
                for record in self._mesh_records
                if record["name"].startswith(f"{name}_")
            ]
            if not matches:
                raise KeyError(f"No surface named {name!r} exists")
            return matches[-1]
        panel = self.active_panel
        for record in reversed(self._mesh_records):
            if record["panel"] == panel:
                return record
        raise RuntimeError(f"No surface exists in panel {panel}")

    def surface(
        self,
        vertices: Any,
        triangles: Any,
        name: str | None = None,
        enable_print: bool = False,
        color: str | None = None,
        opacity: float | None = None,
    ) -> None:
        self._ensure_started()
        if isinstance(vertices, (list, tuple)) and isinstance(triangles, (list, tuple)):
            vertices, triangles = self._combine_meshes(vertices, triangles)
        if self._should_debounce():
            return

        mesh = self._make_polydata(vertices, triangles)
        object_name = (
            f"obj_{self._marker_counter}"
            if name is None
            else f"{name}_{self._marker_counter}"
        )
        self._marker_counter += 1
        record = {
            "name": object_name,
            "mesh": mesh,
            "panel": self.active_panel,
            "actor": None,
            "color": self._normalize_color(color or "lightgray"),
            "opacity": self._opacity if opacity is None else float(opacity),
            "show_edges": self._edge_color is not None,
            "color_limits": None,
        }
        self._mesh_records.append(record)
        self._objects.append(object_name)
        self._draw_mesh(record)
        if enable_print:
            print(f"Surface: {mesh.n_points} verts, {mesh.n_cells} cells")

    def values(
        self, fun: Any, name: str = "", vmin: float | None = None, vmax: float | None = None
    ) -> None:
        self._ensure_started()
        if self._should_debounce():
            return
        values = np.asarray(fun, dtype=float)
        if values.ndim == 2:
            if values.shape[1] == 0:
                raise ValueError("fun must contain at least one value column")
            series = values
        elif values.ndim != 1:
            raise ValueError("fun must have shape (N,) or (N, M)")
        else:
            series = values.reshape(-1, 1)
        record = self._find_mesh(name)
        if series.shape[0] != record["mesh"].n_points:
            raise ValueError(
                f"fun has {series.shape[0]} values but surface has "
                f"{record['mesh'].n_points} vertices"
            )
        record["series"] = series
        values = series[:, min(self._frame_index, series.shape[1] - 1)]
        record["mesh"].point_data["values"] = values
        self._refresh_data_range()
        if vmin is not None or vmax is not None:
            data_limits = self._data_range or (0.0, 1.0)
            self.color_range(
                data_limits[0] if vmin is None else vmin,
                data_limits[1] if vmax is None else vmax,
            )
        elif not self._manual_color_limits and self._data_range is not None:
            self._apply_color_range(self._data_range)
        elif self._color_limits is not None and self._lookup_table is None:
            self._apply_color_range(self._color_limits)
        record["color_limits"] = self._color_limits
        self._remove_actor(record)
        self._draw_mesh(record)
        self._create_scalar_bar()
        self._update_controls()
        if self._picked_node is not None and self._picked_node[0] is record:
            self._update_picked_node_info()
        self.plotter.render()

    def set_frame(self, frame: int) -> None:
        """Select a scalar-data column for every surface that has time frames."""
        frame = int(frame)
        frame_counts = [
            record["series"].shape[1]
            for record in self._mesh_records
            if record.get("series") is not None
        ]
        if not frame_counts:
            raise RuntimeError("No scalar frames are available")
        max_frame = max(frame_counts) - 1
        if frame < 0 or frame > max_frame:
            raise ValueError(f"frame must be between 0 and {max_frame}")
        self._frame_index = frame
        for record in self._mesh_records:
            series = record.get("series")
            if series is None:
                continue
            record["mesh"].point_data["values"] = series[
                :, min(frame, series.shape[1] - 1)
            ]
            self._remove_actor(record)
            self._draw_mesh(record)
        if self.plotter is not None:
            self._update_controls()
            self._update_picked_node_info()
            self.plotter.render()

    @staticmethod
    def _validate_bins(bins: int) -> None:
        if int(bins) != bins or bins < 1:
            raise ValueError("bins must be a positive integer")

    def _refresh_data_range(self) -> None:
        finite_arrays = [
            record["series"][np.isfinite(record["series"])]
            for record in self._mesh_records
            if record.get("series") is not None
            and np.isfinite(record["series"]).any()
        ]
        if not finite_arrays:
            self._data_range = None
            return
        minimum = min(float(values.min()) for values in finite_arrays)
        maximum = max(float(values.max()) for values in finite_arrays)
        minimum = float(np.floor(minimum))
        maximum = float(np.ceil(maximum))
        if minimum == maximum:
            maximum += 1.0
        self._data_range = (minimum, maximum)

    def _create_lookup_table(self) -> None:
        import pyvista as pv

        self._lookup_table = pv.LookupTable(
            cmap=self._color_map,
            n_values=self._bins,
        )
        if self._color_limits is not None:
            self._lookup_table.scalar_range = self._color_limits
        for record in self._mesh_records:
            if record.get("actor") is not None and "values" in record["mesh"].point_data:
                record["actor"].mapper.lookup_table = self._lookup_table

    def _apply_color_range(self, limits: tuple[float, float]) -> None:
        if limits[0] >= limits[1]:
            raise ValueError("minimum color value must be less than maximum")
        self._color_limits = (float(limits[0]), float(limits[1]))
        self.minmax_values = self._color_limits
        if self._lookup_table is None:
            self._create_lookup_table()
        self._lookup_table.scalar_range = self._color_limits
        for record in self._mesh_records:
            if "values" in record["mesh"].point_data:
                record["color_limits"] = self._color_limits
                actor = record.get("actor")
                if actor is not None:
                    actor.mapper.lookup_table = self._lookup_table
                    actor.mapper.scalar_range = self._color_limits
        if self.plotter is not None:
            if self._scalar_bar_added:
                self.plotter.update_scalar_bar_range(
                    self._color_limits,
                    name="Result values",
                )
            self.plotter.render()
        self._sync_display_fields()

    def _create_scalar_bar(self) -> None:
        if self.plotter is None or self._scalar_bar_added or self._lookup_table is None:
            return
        record = next(
            (
                item
                for item in self._mesh_records
                if item.get("actor") is not None and "values" in item["mesh"].point_data
            ),
            None,
        )
        if record is None:
            return
        self._select_panel((1, 1))
        self.plotter.add_scalar_bar(
            title="Result values",
            mapper=record["actor"].mapper,
            n_labels=5,
            n_colors=self._bins,
            vertical=False,
            width=0.5,
            height=0.08,
            position_x=0.25,
            position_y=0.88,
            interactive=False,
            render=False,
        )
        self._scalar_bar_added = True

    def _update_controls(self) -> None:
        if self.plotter is None:
            return
        self._select_panel((1, 1))
        self._updating_controls = True
        try:
            self.plotter.clear_slider_widgets()
            frame_counts = [
                record["series"].shape[1]
                for record in self._mesh_records
                if record.get("series") is not None
            ]
            if frame_counts and max(frame_counts) > 1:
                max_frame = max(frame_counts) - 1
                self._frame_index = min(self._frame_index, max_frame)
                self.plotter.add_slider_widget(
                    self._on_frame_slider,
                    rng=(0, max_frame),
                    value=self._frame_index,
                    title="Frame",
                    pointa=(0.3, 0.16),
                    pointb=(0.7, 0.16),
                )
            self._sync_display_fields()
            self._select_panel(self.active_panel)
        finally:
            self._updating_controls = False

    def _on_frame_slider(self, value: float) -> None:
        if not self._updating_controls:
            self.set_frame(round(value))

    def _apply_display_settings(self, checked: bool = False) -> None:
        del checked
        if (
            self._range_min_field is None
            or self._range_max_field is None
            or self._bins_field is None
            or self._opacity_field is None
        ):
            return
        try:
            minimum = float(self._range_min_field.text())
            maximum = float(self._range_max_field.text())
            bins = int(self._bins_field.text())
            opacity = float(self._opacity_field.text())
            self._validate_bins(bins)
            if not np.isfinite((minimum, maximum)).all():
                raise ValueError("Color range values must be finite.")
            if minimum >= maximum:
                raise ValueError("Color minimum must be less than color maximum.")
            if not np.isfinite(opacity) or not 0 <= opacity <= 1:
                raise ValueError("Opacity must be between 0 and 1.")
        except ValueError as exc:
            self._controls_status.setText(str(exc))
            return

        self.color_range(minimum, maximum)
        self.set_bins(bins)
        self.transparency(opacity)
        self._controls_status.setText("Display settings applied.")

    def set_transparency_mode(self, mode: str) -> None:
        """Choose standard alpha blending or depth peeling for translucent surfaces."""
        normalized = mode.strip().lower().replace("-", "_").replace(" ", "_")
        if normalized in {"standard", "alpha", "standard_alpha_blending"}:
            normalized = "standard"
        elif normalized in {"depth_peeling", "depthpeeling"}:
            normalized = "depth_peeling"
        else:
            raise ValueError("mode must be 'standard' or 'depth_peeling'")

        if normalized == self._transparency_mode and self.plotter is not None:
            return
        previous_mode = self._transparency_mode
        self._transparency_mode = normalized
        if self.plotter is not None:
            try:
                self._apply_transparency_mode()
            except RuntimeError:
                self._transparency_mode = previous_mode
                self._sync_display_fields()
                raise
        self._sync_display_fields()
        if self.plotter is not None:
            self.plotter.render()

    def _apply_transparency_mode(self) -> None:
        if self.plotter is None:
            return
        if self._transparency_mode == "standard":
            for renderer in self.plotter.renderers:
                renderer.disable_depth_peeling()
            return
        unsupported_renderers = [
            renderer
            for renderer in self.plotter.renderers
            if not renderer.enable_depth_peeling(
                number_of_peels=8,
                occlusion_ratio=0.0,
            )
        ]
        if unsupported_renderers:
            for renderer in self.plotter.renderers:
                renderer.disable_depth_peeling()
            raise RuntimeError(
                "Depth peeling is not supported by this graphics driver. "
                "Use Standard alpha blending instead."
            )

    def _remove_actor(self, record: dict[str, Any]) -> None:
        actor = record.get("actor")
        if actor is not None:
            self.plotter.remove_actor(actor)
            record["actor"] = None

    def color_range(self, vmin: float, vmax: float) -> None:
        limits = (float(vmin), float(vmax))
        if limits[0] >= limits[1]:
            raise ValueError("vmin must be less than vmax")
        self._manual_color_limits = True
        self._apply_color_range(limits)
        self._create_scalar_bar()
        self._update_controls()

    def gradient_bins(self, bins: int = 15) -> None:
        """Set the number of discrete colors used in the shared scalar map."""
        self._validate_bins(bins)
        self._bins = int(bins)
        self._replace_lookup_table()
        self._refresh_scalar_meshes()

    def set_bins(self, bins: int = 15) -> None:
        """Set the number of discrete colors used in the shared scalar map."""
        self.gradient_bins(bins)
        self._sync_display_fields()

    def _replace_lookup_table(self) -> None:
        if self._lookup_table is None:
            return
        self._remove_scalar_bar()
        self._create_lookup_table()
        self._create_scalar_bar()

    def _remove_scalar_bar(self) -> None:
        if self.plotter is not None and self._scalar_bar_added:
            self.plotter.remove_scalar_bar(title="Result values", render=False)
        self._scalar_bar_added = False

    def _refresh_scalar_meshes(self) -> None:
        if self.plotter is None:
            return
        self._remove_scalar_bar()
        for record in self._mesh_records:
            if "values" in record["mesh"].point_data:
                self._remove_actor(record)
                self._draw_mesh(record)
        self._create_scalar_bar()
        self._update_controls()
        self.plotter.render()

    def set_active_panel(self, horizontal: int = 1, vertical: int = 1) -> None:
        horizontal, vertical = int(horizontal), int(vertical)
        self.active_panel = (
            min(max(horizontal, 1), self._panel_count[0]),
            min(max(vertical, 1), self._panel_count[1]),
        )
        self._select_panel(self.active_panel)

    def set_panels_number(self, horizontal: int = 1, vertical: int = 1) -> None:
        horizontal, vertical = int(horizontal), int(vertical)
        if horizontal < 1 or vertical < 1:
            raise ValueError("panel counts must be positive integers")
        if (horizontal, vertical) == self._panel_count:
            return
        self._panel_count = (horizontal, vertical)
        self.panels = self._panel_count
        self.active_panel = (
            min(self.active_panel[0], horizontal),
            min(self.active_panel[1], vertical),
        )
        for record in self._mesh_records:
            record["panel"] = (
                min(record["panel"][0], horizontal),
                min(record["panel"][1], vertical),
            )
        for record in self._text_records:
            record["panel"] = (
                min(record["panel"][0], horizontal),
                min(record["panel"][1], vertical),
            )
        self._refresh_plotter()

    def text(self, text: str, pos: tuple[float, float] = (0, 0)) -> None:
        if any(value < 0 or value > 1 for value in pos):
            raise ValueError("pos must be in (0,1) range for both x and y")
        self._ensure_started()
        record = {
            "text": str(text),
            "position": tuple(pos),
            "color": self._text_color,
            "font_size": 16,
            "panel": self.active_panel,
            "actor": None,
        }
        self._text_records.append(record)
        self._draw_text(record)

    def property_on_mouse_click(self, property_name: str) -> None:
        supported = {"coor", "tri", "ver", "fun", "lm", "cross"}
        if property_name not in supported:
            raise ValueError(f"property_name must be one of {sorted(supported)}")
        self._pick_property = property_name
        self._ensure_started()
        self.plotter.enable_point_picking(
            callback=self._on_pick,
            show_message=False,
            use_picker=True,
        )

    def _on_pick(self, point: Any, picker: Any = None) -> None:
        if point is None:
            return
        point = np.asarray(point, dtype=float)
        dataset = picker.GetDataSet() if picker is not None else None
        record = next(
            (
                item
                for item in self._mesh_records
                if dataset is not None and dataset == item["mesh"]
            ),
            None,
        )
        if record is None:
            record = self._find_mesh()
        mesh = record["mesh"]
        point_id = int(picker.GetPointId()) if picker is not None else -1
        if point_id < 0 or point_id >= mesh.n_points:
            point_id = int(mesh.find_closest_point(point))
        self._picked_node = (record, point_id)
        self._update_picked_node_info()

    def _reset_picked_node_info(self) -> None:
        self._picked_node = None
        for label in (
            self._picked_node_number,
            self._picked_node_coordinates,
            self._picked_node_value,
        ):
            if label is not None:
                label.setText("—")

    def _update_picked_node_info(self) -> None:
        if self._picked_node is None:
            return
        record, point_id = self._picked_node
        if (
            all(item is not record for item in self._mesh_records)
            or point_id < 0
            or point_id >= record["mesh"].n_points
        ):
            self._reset_picked_node_info()
            return

        coordinates = record["mesh"].points[point_id]
        series = record.get("series")
        if series is None:
            value_text = "N/A"
        else:
            frame = min(self._frame_index, series.shape[1] - 1)
            value = float(series[point_id, frame])
            value_text = f"{value:.8g}" if np.isfinite(value) else "N/A"

        if self._picked_node_number is not None:
            self._picked_node_number.setText(str(point_id + 1))
            self._picked_node_coordinates.setText(
                ", ".join(f"{coordinate:.8g}" for coordinate in coordinates)
            )
            self._picked_node_value.setText(value_text)

    def reset(self) -> None:
        self.clear()
        if self.plotter is not None:
            self.plotter.reset_camera()
            self.plotter.render()

    def clear(self, enable_print: bool = False) -> None:
        del enable_print
        if self.plotter is not None:
            self.plotter.clear_slider_widgets()
            self._remove_scalar_bar()
            self.plotter.clear_actors()
            self.plotter.add_axes()
        self._mesh_records.clear()
        self._text_records.clear()
        self._objects.clear()
        self._marker_counter = 0
        self._frame_index = 0
        self._data_range = None
        self._lookup_table = None
        self._reset_picked_node_info()

    def close(self) -> None:
        if self.plotter is not None:
            self._dispose_controls_dock()
            self.plotter.close()
            self.plotter = None

    def cmd(self, cmd: str, enable_print: bool = False) -> None:
        self._ensure_started()
        if enable_print:
            print("Command:", cmd)
        parts = cmd.strip().split(maxsplit=3)
        if not parts:
            return
        command = parts[0].lower()
        args = parts[1:]
        if command == "reset":
            self.reset()
        elif command in {"delete", "del"}:
            if not args:
                raise ValueError("delete requires an object name")
            matches = [item for item in self._mesh_records if item["name"] == args[0]]
            for record in matches:
                self._remove_actor(record)
                self._mesh_records.remove(record)
                if record["name"] in self._objects:
                    self._objects.remove(record["name"])
                if self._picked_node is not None and self._picked_node[0] is record:
                    self._reset_picked_node_info()
        elif command == "trans":
            if not args:
                raise ValueError("trans requires an opacity")
            opacity = float(args[0])
            if not 0 <= opacity <= 1:
                raise ValueError("opacity must be between 0 and 1")
            self._opacity = opacity
            for record in self._mesh_records:
                record["opacity"] = opacity
                self._remove_actor(record)
                self._draw_mesh(record)
            self._sync_display_fields()
        elif command in {"transmode", "rendering"}:
            if not args:
                raise ValueError("transmode requires standard or depth_peeling")
            self.set_transparency_mode(args[0])
        elif command == "edge":
            if not args:
                raise ValueError("edge requires a color")
            self._edge_color = (
                None
                if args[0].lower() in {"off", "n", "none"}
                else self._normalize_color(args[0])
            )
            for record in self._mesh_records:
                record["show_edges"] = self._edge_color is not None
                self._remove_actor(record)
                self._draw_mesh(record)
        elif command in {"bgdcolor", "background"}:
            if not args:
                raise ValueError("bgdcolor requires a color")
            self._background = self._normalize_color(args[0])
            self.plotter.set_background(self._background)
        elif command == "color":
            if not args:
                raise ValueError("color requires a color")
            record = self._find_mesh()
            record["color"] = self._normalize_color(args[0])
            self._remove_actor(record)
            self._draw_mesh(record)
        elif command in {"funscale", "clim"}:
            if len(args) != 2:
                raise ValueError("funscale requires a minimum and maximum")
            self.color_range(float(args[0]), float(args[1]))
        elif command == "step":
            if not args:
                raise ValueError("step requires a gradient step")
            step = float(args[0])
            if step <= 0:
                raise ValueError("gradient step must be positive")
            if self._color_limits is None:
                raise RuntimeError("Set or display a scalar range before using step")
            bins = int(np.ceil((self._color_limits[1] - self._color_limits[0]) / step))
            self.set_bins(max(1, bins))
        elif command == "textcolor":
            if not args:
                raise ValueError("textcolor requires a color")
            self._text_color = args[0]
        elif command == "text":
            if len(parts) < 4:
                raise ValueError("text requires x, y, and text")
            self.text(parts[3], (float(args[0]), float(args[1])))
        elif command == "panel":
            if len(args) != 2:
                raise ValueError("panel requires horizontal and vertical indices")
            self.set_active_panel(int(args[0]), int(args[1]))
        elif command == "horizontal":
            if not args:
                raise ValueError("horizontal requires a panel count")
            self.set_panels_number(int(args[0]), self._panel_count[1])
        elif command == "vertical":
            if not args:
                raise ValueError("vertical requires a panel count")
            self.set_panels_number(self._panel_count[0], int(args[0]))
        elif command == "mouse":
            if not args:
                raise ValueError("mouse requires a property name")
            self.property_on_mouse_click(args[0])
        elif command == "angle":
            if len(args) != 3:
                raise ValueError("angle requires x, y, and z angles")
            camera = self.plotter.camera
            camera.Azimuth(float(args[0]))
            camera.Elevation(float(args[1]))
            camera.Roll(float(args[2]))
            self.plotter.render()
        else:
            raise ValueError(f"Unsupported QTriplot command: {cmd!r}")

    @staticmethod
    def _normalize_color(color: str) -> str:
        return {
            "r": "red",
            "g": "green",
            "b": "blue",
            "c": "cyan",
            "m": "magenta",
            "y": "yellow",
            "k": "black",
            "w": "white",
        }.get(color.lower(), color)


__all__ = ["Result3D"]
