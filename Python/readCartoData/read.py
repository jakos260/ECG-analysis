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
from Carto3Reader.CartoMap import CartoMap

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

v = data['vertices']
t = data['triangles']
points = data['points']
colors_mesh = data['colors_mesh']
color_names = data['color_names']


# Drop columns with zero standard deviation
if colors_mesh is not None and color_names:
    # stds = np.nanstd(colors_mesh, axis=0)
    valid_indices = np.min(colors_mesh, axis=0) != np.max(colors_mesh, axis=0)
    colors_mesh = colors_mesh[:, valid_indices]
    color_names = [name for name, valid in zip(color_names, valid_indices) if valid]

# Plot histograms for all color types on a single chart with transparency
# if colors_mesh is not None and color_names:
#     fig = go.Figure()
    
#     for i, name in enumerate(color_names):
#         values = colors_mesh[:, i]
#         valid_values = values[~np.isnan(values)]
#         if len(valid_values) > 0:
#             fig.add_trace(go.Histogram(x=valid_values, nbinsx=30, name=name, opacity=0.5))
    
#     fig.update_layout(
#         title_text="Histograms of Mesh Color Values for All Types",
#         barmode='overlay'
#     )
#     fig.show()

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
# unipolar_mesh = np.where(unipolar_mesh < 20, unipolar_mesh, 10)  # drop values above 20mV (likely noise) to 0 for better visualization

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
q.text(f"Carto3Data {measurement_name} - {map_type}", pos=(0.25, 0.95))
# q.background_color("white")
q.markers([(p.X, p.Y, p.Z) for p in points.itertuples()], color='red', r=1)

input("Press Enter to close QTriplot...")
q.close()

# fig = go.Figure()
# if v is not None and t is not None:
#     fig.add_trace(go.Mesh3d(
#         x=v[:, 0], y=v[:, 1], z=v[:, 2],
#         i=t[:, 0], j=t[:, 1], k=t[:, 2],
#         color='lightgray',
#         opacity=0.3, # Półprzezroczystość, żeby było widać wewnątrz
#         name='Anatomia'
#     ))

# fig.add_trace(go.Scatter3d(
#     x=df['X'], y=df['Y'], z=df['Z'],
#     mode='markers',
#     marker=dict(
#         size=5,
#         color=df['Bipolar'],
#         colorscale='Jet',
#         colorbar=dict(title="Woltaż [mV]"),
#         cmin=0, cmax=1.5, # Zazwyczaj w Carto granica zdrowej tkanki to > 1.5mV
#         opacity=0.9
#     ),
#     text="ID: " + df['ID'].astype(str) + "<br>V: " + df['Bipolar'].astype(str) + " mV",
#     name='Punkty Mapowania',
#     hoverinfo='text'
# ))

# fig.update_layout(
#     scene=dict(bgcolor='black', xaxis_visible=False, yaxis_visible=False, zaxis_visible=False),
#     title="Model 3D - Lewa Komora (Carto 3)"
# )
# fig.show()
