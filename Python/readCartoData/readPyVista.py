import os
import sys
import numpy as np
import pyvista as pv
from pathlib import Path
from dotenv import load_dotenv

# Konfiguracja ścieżek
root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))                         
from carto3reader.CartoMap import CartoMap

load_dotenv()

data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
patient_name = 'Patient 2025_10_20'
sample_data_path = os.path.join(data_path, 'other', 'Carto3Data', patient_name, 'VT1', 'Export_VT1-02_02_2026-12-24-27')

measurement_names = [
    '1-VT' 
    # '2-AT LV', 
    # '2-1-ReAT LV', 
    # '1-1-1-PATTERN-STYMU', 
    # '1-1-PATTERN BEZ STYMU'
]

default_map_type = 'Unipolar'
plotter = pv.Plotter(title=f"Carto3Data - {patient_name} {default_map_type}")

for name in measurement_names:
    print(f"Loading model: {name}...")
    try:
        carto_map = CartoMap(sample_data_path, name)
        data = carto_map.load_mesh()

        # 1. Vertices validation and transformation
        v = np.array(data['vertices'])
        
        if v.shape[0] == 3 and v.shape[1] > 3:
            print(f"  -> [{name}] Transposing vertices from {v.shape} to (N, 3)")
            v = v.T
            
        v = np.ascontiguousarray(v, dtype=np.float64)

        # 2. Triangles validation and transformation
        t = np.array(data['triangles'])
        
        if t.shape[0] == 3 and t.shape[1] > 3:
            print(f"  -> [{name}] Transposing triangles from {t.shape} to (M, 3)")
            t = t.T

        # --- FIX: Filter out invalid triangles with negative indices ---
        initial_t_count = t.shape[0]
        # Keep only rows where ALL indices are >= 0
        valid_mask = (t >= 0).all(axis=1)
        t = t[valid_mask]
        filtered_count = initial_t_count - t.shape[0]
        
        if filtered_count > 0:
            print(f"  -> [{name}] Removed {filtered_count} invalid triangles (negative indices found).")

        # 3. Dynamic Indexing correction (0-based vs 1-based)
        if t.size > 0:
            t_min = t.min()
            print(f"  -> [{name}] Valid minimum triangle index is {t_min}")
            
            if t_min == 1:
                print(f"  -> [{name}] Converting 1-based indexing to 0-based")
                t = t - 1
        else:
            print(f"  -> [{name}] ERROR: No valid triangles left after filtering!")
            continue # Skip rendering this model

        t = np.ascontiguousarray(t, dtype=np.int32)

        # 4. PyVista faces array construction
        faces = np.insert(t, 0, 3, axis=1).ravel()

        # 5. Create Mesh object
        mesh_pv = pv.PolyData(v, faces)

        # 6. Mesh coloring
        colors_mesh = data['colors_mesh']
        color_names = data['color_names']
        
        map_type = default_map_type if default_map_type in color_names else (color_names[0] if color_names else None)

        if map_type is not None and colors_mesh is not None:
            color_idx = color_names.index(map_type)
            mesh_values = colors_mesh[:, color_idx]
            
            mesh_pv.point_data[map_type] = mesh_values
            
            plotter.add_mesh(mesh_pv, scalars=map_type, cmap='jet', opacity=0.8, 
                             show_scalar_bar=True, label=f"{name}", show_edges=True)
        else:
            print(f"  -> Warning: No valid color data for {name}, displaying with solid color.")
            plotter.add_mesh(mesh_pv, color='lightgray', opacity=0.6, label=f"{name} (no colors)", show_edges=True)

        # 7. Electrodes rendering
        cs = carto_map.load_electrodes('CS')
        if cs:
            cs_points = np.array([(cs[k][0][1], cs[k][0][2], cs[k][0][3]) for k in cs.keys()])
            if cs_points.size > 0:
                plotter.add_points(pv.PolyData(cs_points), color='red', point_size=10, render_points_as_spheres=True)

    except Exception as e:
        print(f"Error: Failed to load or render {name}. Details: {e}")

plotter.add_legend()
plotter.add_axes()
plotter.set_background('white')

plotter.show()