import os, sys
import numpy as np
import pyvista as pv
from matplotlib.colors import ListedColormap

from pathlib import Path
from dotenv import load_dotenv

root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))                         # add helpers
from qtripy.readECGsim import read_ecg_sim
from helpers.split_aha_segments import get_heart_model_segments

load_dotenv()
data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
sample_data_path = os.path.join(data_path, "ECGsim_data", "normal_young_male")

heart_data = read_ecg_sim(sample_data_path)
ventricles_ver, ventriclec_tri = heart_data['VENTR']['geom']['VER'], heart_data['VENTR']['geom']['ITRI']

apex_idx = 578              # Wierzchołek koniuszka
base_idx = 513              # Wierzchołek gdzieś na środku podstawy zastawek
septum_bottom_idx = 564     # Wierzchołek na dole przegrody
septum_top_idx = 995        # Wierzchołek na górze przegrody 

apex_pt = ventricles_ver[apex_idx]
base_pt = ventricles_ver[base_idx]
septum_bottom_pt = ventricles_ver[septum_bottom_idx]
septum_top_pt = ventricles_ver[septum_top_idx]

aha_labels = get_heart_model_segments(
    ventricles_ver,
    apex_pt,
    base_pt,
    septum_top_pt,
    septum_bottom_pt,
    angle_offset_deg=0,
    extra_septum=True
    )

unique_labels_number = np.unique(aha_labels).shape[0]
print(f"labeled {aha_labels.shape[0]} vertices with AHA {unique_labels_number} unique segments")


#visualization
ventriclec_tri = ventriclec_tri.astype(int)
faces = np.column_stack((np.full(len(ventriclec_tri), 3), ventriclec_tri)).flatten()

mesh = pv.PolyData(ventricles_ver, faces)
tri_labels = aha_labels[ventriclec_tri[:, 0]]
mesh.cell_data["Segments"] = tri_labels

colors = [
    "white",      
    "purple",     
    "grey",       
    "black",      
    "red",        
    "orange",  
    "#CC1933", 
    "green",      
    "yellow",     
    "cyan",       
    "purple",     
    "#FF8080",   
    "darkgrey",   
    "lightblue",  
    "lightyellow",
    "yellow",     
    "black",     
    "#FFCC99",
    "blue",
    "lime",  
]
my_cmap = ListedColormap(colors[:unique_labels_number+1])

plotter = pv.Plotter()
plotter.add_mesh(
    mesh,
    scalars="Segments",
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

file_name = "ventricles_labels_of_20.segments"
np.savetxt(os.path.join(sample_data_path, "model", file_name), aha_labels, fmt="%d")

with open(os.path.join(sample_data_path, "model", "ventricles_labels_points.segments"), "w") as f:
    lines = [
        "'file_name','apex','base','septum_top','septum_bottom'\n",
        f"{file_name},{apex_idx},{base_idx},{septum_top_idx},{septum_bottom_idx}\n"]
    f.writelines(lines)