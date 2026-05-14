import os, sys
import numpy as np
import pyvista as pv
from matplotlib.colors import ListedColormap

from pathlib import Path
from dotenv import load_dotenv

root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))                         # add helpers
from qtripy.readECGsim import read_ecg_sim
from helpers.split_aha_segments import get_17_aha_segments

load_dotenv()
data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
sample_data_path = os.path.join(data_path, "ECGsim_data", "normal_young_male")

heart_data = read_ecg_sim(sample_data_path)
ventricles_ver, ventriclec_tri = heart_data['VENTR']['geom']['VER'], heart_data['VENTR']['geom']['ITRI']

apex_idx = 558                 # Indeks wierzchołka na koniuszku (np. linia 105 w pliku)
septum_idx = 249            # Wierzchołek gdzieś na środku podstawy zastawek

basal_rim_indices = [532, 513, 1019, 1011, 505] # Obwód zastawki mitralnej
rim_points = ventricles_ver[basal_rim_indices]

base_pt = np.mean(rim_points, axis=0)
apex_pt = ventricles_ver[apex_idx]
septum_pt = ventricles_ver[septum_idx]


aha_labels = get_17_aha_segments(
    ventricles_ver,
    apex_pt,
    base_pt,
    septum_pt,
    angle_offset_deg=0,
    extra_septum=True,
    septum_wedge_deg=40
    )
unique_labels_number = np.unique(aha_labels).shape[0]
print(f"labeled {aha_labels.shape[0]} vertices with AHA {unique_labels_number} unique segments")


#visualization
ventriclec_tri = ventriclec_tri.astype(int)
faces = np.column_stack((np.full(len(ventriclec_tri), 3), ventriclec_tri)).flatten()

mesh = pv.PolyData(ventricles_ver, faces)
tri_labels = aha_labels[ventriclec_tri[:, 0]]
mesh.cell_data["AHA_Segment"] = tri_labels

colors = [
    "white",           # 1
    "purple",          # 11
    "grey",            # 3
    "black",           # 5
    "red",             # 6
    "lightgrey",           # 7
    "#CC1933",        # 17 (heart - bordowy)
    "green",       # 8
    "yellow",          # 9
    "cyan",            # 10
    "purple",          # 11
    "#FF8080",       # 12 (lightred)
    "darkgrey",        # 4
    # "lightgreen",      # 13
    "lightblue",       # 14
    "lightyellow",     # 15
    "yellow",       # 14
    "black",     # 15
    "#FFCC99",       # 16 (flesh - cielisty)
    "blue",
    "lime"              # 2
]
my_cmap = ListedColormap(colors)

plotter = pv.Plotter()
plotter.add_mesh(
    mesh,
    scalars="AHA_Segment",
    cmap=my_cmap,
    show_scalar_bar=True,
    show_edges=True,
    scalar_bar_args={
        "title": "Segment AHA",
        "fmt": "%.0f",
        "n_labels": unique_labels_number
    }
)

def pick_callback(point):
    vertex_id = mesh.find_closest_point(point)
    print(f"Zaznaczono wierzchołek ID: {vertex_id:5} | Współrzędne: {np.round(point, 2)}")

plotter.enable_point_picking(
    callback=pick_callback,
    left_clicking=True,
    show_point=True, 
    point_size=10
)

plotter.set_background("#9999CC") 
plotter.show()
