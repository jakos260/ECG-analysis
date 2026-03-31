import os, sys
import numpy as np
import xml.etree.ElementTree as ET
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

from pathlib import Path
from dotenv import load_dotenv

root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))                         # add helpers
from qtripy.qtripy import QTripy
from carto3reader.CartoMap import CartoMap

load_dotenv()

data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
sample_data_path = os.path.join(data_path, 'Carto3Data', 'Patient 2025_10_20', 'VT1', 'Export_VT1-02_02_2026-12-24-27')

measurement_name = '1-VT' # 114 points
# measurement_name = '2-AT LV'
# measurement_name = '2-1-ReAT LV'
# measurement_name = '1-1-1-PATTERN-STYMU'
# measurement_name = '1-1-PATTERN BEZ STYMU'

carto_map = CartoMap(sample_data_path, measurement_name)
data = carto_map.load_mesh()

cs = carto_map.load_electrodes('CS')
mc = carto_map.load_electrodes('MCC_DX')
qd = carto_map.load_electrodes('QUAD_A')

v = data['vertices']
t = data['triangles']
points = data['points']
colors_mesh = data['colors_mesh']
color_names = data['color_names']


# Drop columns with zero standard deviation
if colors_mesh is not None and color_names:
    valid_indices = np.min(colors_mesh, axis=0) != np.max(colors_mesh, axis=0)
    colors_mesh = colors_mesh[:, valid_indices]
    color_names = [name for name, valid in zip(color_names, valid_indices) if valid]

def choose_color_type(color_names, default='Unipolar'):
    if not color_names:
        return None

    print("Available colorings:")
    for i, name in enumerate(color_names, start=1):
        marker = " (default)" if name == default else ""
        print(f"  {i}. {name}{marker}")

    default_index = color_names.index(default) + 1 if default in color_names else 1
    prompt = f"Choose color index [1-{len(color_names)}] (default {default_index}): "

    try:
        import msvcrt
        print("Press number key or Enter to select default.")
        while True:
            ch = msvcrt.getwch()
            if ch in ('\r', '\n'):
                return default if default in color_names else color_names[0]
            if ch.isdigit():
                idx = int(ch)
                if 1 <= idx <= len(color_names):
                    print(ch)
                    return color_names[idx - 1]
                else:
                    print(f"Invalid choice {ch}, expected 1-{len(color_names)}")
    except Exception:
        # Fallback for non-Windows or if msvcrt can't be used
        choice = input(prompt).strip()
        if not choice:
            return default if default in color_names else color_names[0]
        if choice.isdigit():
            idx = int(choice)
            if 1 <= idx <= len(color_names):
                return color_names[idx - 1]
        print("Invalid selection, using default.")
        return default if default in color_names else color_names[0]


map_type = choose_color_type(color_names, default='Unipolar')
mesh = colors_mesh[:, color_names.index(map_type)] if (map_type is not None and map_type in color_names) else None

q = QTripy()
q.begin()
q.reset()
# q.set_panels_number(2,1)
q.surface(v, t, set_max_distance=True)
q.transparency(0.5)
q.values(mesh)
q.gradient_bins(10)

# q.set_active_panel(2,1)
# q.surface(v, t)
# q.transparency(0.3)

q.property_on_mouse_click('coor')
q.text(f"Carto3Data {measurement_name} - {map_type}", pos=(0.15, 0.95))
# q.background_color("white")

cs_points = [(cs[lead_key][0][1], cs[lead_key][0][2], cs[lead_key][0][3]) for lead_key in cs.keys()]
mc_points = [(mc[lead_key][0][1], mc[lead_key][0][2], mc[lead_key][0][3]) for lead_key in mc.keys()]
qd_points = [(qd[lead_key][0][1], qd[lead_key][0][2], qd[lead_key][0][3]) for lead_key in qd.keys()]

q.markers(cs_points, color='red', r=1)
q.markers(mc_points, color='blue', r=1)
q.markers(qd_points, color='green', r=1)


# carto_map.ecg_reader.plot()  # Plot ECG data if available

input("Press Enter to close QTriplot...")
q.close()
