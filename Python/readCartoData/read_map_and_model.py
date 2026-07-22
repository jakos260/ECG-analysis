import os
import sys
import numpy as np
import pyvista as pv
from pathlib import Path
from dotenv import load_dotenv

# Path configuration
root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))

from carto3reader.CartoMap import CartoMap
from helpers.reader import loadtri

load_dotenv()


def make_polydata_mesh(vertices: np.ndarray, triangles: np.ndarray):
    if vertices.size == 0:
        return None

    faces = np.hstack((np.full((triangles.shape[0], 1), 3), triangles)).ravel()
    return pv.PolyData(vertices, faces)


model = dict()

data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
model_data_path = os.path.join(data_path, 'other', 'ModelDataBase', 'normalM40')
carto_data_path = os.path.join(data_path, 'other', 'Carto3Data', 'Patient 2025_10_20', 'VT1', 'Export_VT1-02_02_2026-12-24-27')

measurement_names = ['1-VT']
default_map_type = 'Unipolar'

plotter = pv.Plotter(title="Carto3Data Models Visualization - PyVista")

for name in measurement_names:
    print(f"Loading model: {name}...")
    try:
        model["lcav"] = loadtri(os.path.join(model_data_path, 'normalM40_lcav.tri'))
        model["rcav"] = loadtri(os.path.join(model_data_path, 'normalM40_rcav.tri'))

        for k, (vertices, triangles) in model.items():
            cavity_mesh = make_polydata_mesh(vertices, triangles)
            if cavity_mesh is None:
                continue

            plotter.add_mesh(
                cavity_mesh,
                color='cyan' if k == 'lcav' else 'magenta',
                opacity=0.55,
                show_edges=True,
                label=f"{k}"
            )

        carto_map = CartoMap(carto_data_path, name)
        data = carto_map.load_mesh()

        v = data['vertices']
        t = np.array(data['triangles']) # Conversion indexing from Fortran 1-based to Python 0-based

        colors_mesh = data['colors_mesh']
        color_names = data['color_names']

        # 1. Triangle conversion for PyVista
        # Format: [3, v1, v2, v3, 3, v1, v2, v3...]
        faces = np.hstack((np.full((t.shape[0], 1), 3), t)).ravel()

        mesh_pv = pv.PolyData(v, faces)
        map_type = default_map_type if default_map_type in color_names else (color_names[0] if color_names else None)

        if map_type is not None and colors_mesh is not None:
            color_idx = color_names.index(map_type)
            mesh_values = colors_mesh[:, color_idx]
            mesh_pv.point_data[map_type] = mesh_values
            
            plotter.add_mesh(
                mesh_pv,
                scalars=map_type,
                cmap='jet',
                opacity=0.5,
                show_scalar_bar=True,                
                show_edges=True,
                label=f"{name}"
                )
        else:
            print(f"Warning: No valid color data for {name}, displaying with solid color.")
            plotter.add_mesh(mesh_pv, color='lightgray', opacity=0.6, label=f"{name} (no colors)", show_edges=True)

        # Electrodes rendering block        
        cs = carto_map.load_electrodes('CS')
        mc = carto_map.load_electrodes('MCC_DX')
        qd = carto_map.load_electrodes('QUAD_A')
        if cs:
            plotter.add_points(pv.PolyData(np.array([cs[k][0][1:4] for k in cs.keys()])), color='blue', point_size=10, render_points_as_spheres=True, label='CS Electrodes')
        if mc:
            plotter.add_points(pv.PolyData(np.array([mc[k][0][1:4] for k in mc.keys()])), color='red', point_size=10, render_points_as_spheres=True, label='MCC_DX Electrodes')
        if qd:
            plotter.add_points(pv.PolyData(np.array([qd[k][0][1:4] for k in qd.keys()])), color='green', point_size=10, render_points_as_spheres=True, label='QUAD_A Electrodes')

    except Exception as e:
        print(f"Error: Failed to load or render {name}. Details: {e}")

plotter.add_legend()
plotter.add_axes()
plotter.set_background('white')

plotter.show()