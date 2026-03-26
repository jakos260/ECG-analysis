import os, sys
import xml.etree.ElementTree as ET
import plotly.graph_objects as go
import plotly.express as px

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
faces_type = 'Unipolar'  # 'Unipolar' or 'Bipolar'
carto_map = CartoMap(os.path.join(sample_data_path, f'{measurement_name}_Points_Export.xml'))
carto_map.load_mesh(os.path.join(sample_data_path, f'{measurement_name}.mesh'))

faces = carto_map.to_dataframe()
data = carto_map.get_map_data()
v = data['vertices']
t = data['triangles']
df = data['points_df']
# print(df.head())



faces[faces_type] = faces[faces_type].mask(faces[faces_type] > 20, 0) # drop values above 20mV (likely noise) to 0 for better visualization
# fig = px.histogram(faces[faces_type], nbins=30, title=faces_type)
# fig.show()

q = QTripy()
q.begin()
q.reset()
q.surface(v, t)
q.transparency(0.3)
q.values(faces[faces_type])
q.gradient_bins(10)
q.property_on_mouse_click('coor')
q.text(f"ECGsim {measurement_name} - {faces_type}", pos=(0.25, 0.95))
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
