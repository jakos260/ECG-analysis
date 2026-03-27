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
from CartoMap import CartoMap

load_dotenv()

data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
sample_data_path = os.path.join(data_path, 'Carto3Data', 'Patient 2025_10_20', 'VT1', 'Export_VT1-02_02_2026-12-24-27')

measurement_name = '2-AT LV'
carto_map = CartoMap(os.path.join(sample_data_path, f'{measurement_name}_Points_Export.xml'))

data = carto_map.load_mesh(os.path.join(sample_data_path, f'{measurement_name}.mesh'))

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
    print(f"Kept {len(color_names)} color types with non-zero std: {color_names}")

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


lat_mesh = colors_mesh[:, color_names.index('LAT')] if 'LAT' in color_names else None
bipolar_mesh = colors_mesh[:, color_names.index('Bipolar')] if 'Bipolar' in color_names else None
unipolar_mesh = colors_mesh[:, color_names.index('Unipolar')] if 'Unipolar' in color_names else None
# unipolar_mesh = np.where(unipolar_mesh < 20, unipolar_mesh, 10)  # drop values above 20mV (likely noise) to 0 for better visualization

q = QTripy()
q.begin()
q.reset()
# q.set_panels_number(2,1)
q.surface(v, t)
q.transparency(0.3)
# q.values(unipolar_mesh)
# q.values(bipolar_mesh)
q.values(unipolar_mesh)
q.gradient_bins(10)

# q.set_active_panel(2,1)
# q.surface(v, t)
# q.transparency(0.3)

q.property_on_mouse_click('coor')
q.text(f"ECGsim {measurement_name} - {'unipolar'}", pos=(0.25, 0.95))
# q.background_color("white")

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
