import os, sys
import numpy as np
from pathlib import Path
from dotenv import load_dotenv

root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))                         # add helpers
from qtripy.qtripy import QTripy
from qtripy.readECGsim import read_ecg_sim
from helpers.split_aha_segments import get_17_aha_segments

load_dotenv()
data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
sample_data_path = os.path.join(data_path, "ECGsim_data", "normal_young_male")

heart_data = read_ecg_sim(sample_data_path)
ventricles_ver, ventriclec_tri = heart_data['VENTR']['geom']['VER'], heart_data['VENTR']['geom']['ITRI']

apex_idx = 1497                 # Indeks wierzchołka na koniuszku (np. linia 105 w pliku)
base_center_idx = 76            # Wierzchołek gdzieś na środku podstawy zastawek

basal_rim_indices = [1074, 961, 533, 514, 420] # Obwód zastawki mitralnej
rim_points = ventricles_ver[basal_rim_indices]
septum_pt = np.mean(rim_points, axis=0)

apex_pt = ventricles_ver[apex_idx]
base_pt = ventricles_ver[base_center_idx]

# Alternatywnie możesz podać współrzędne z palca, np:
# base_pt = np.array([47.3, 50.0, -10.0])

aha_labels = get_17_aha_segments(ventricles_ver, apex_pt, base_pt, septum_pt, angle_offset_deg=0)

# output_file = os.path.join(sample_data_path, "model/ventricle_aha_labels.txt")
# np.savetxt(output_file, aha_labels, fmt='%d')
        
colors = [
    "white", "lightgrey", "grey", "darkgrey", "black", "red", "green",
    "blue", "yellow", "cyan", "purple", "lightred", "lightgreen",
    "lightblue", "lightyellow", "flesh", "heart", "slideblue"
]

q = QTripy()
q.begin()
q.reset()


tri_labels = aha_labels[ventriclec_tri[:, 0].astype(int)]

for seg_id in range(1, 18):
    mask = (tri_labels == seg_id)
    segment_triangles = ventriclec_tri[mask]
    
    if len(segment_triangles) > 0:
        current_color = colors[seg_id - 1]
        q.surface(ventricles_ver, segment_triangles, color=current_color)

q.property_on_mouse_click('ver')
q.text("ECGsim Ventricle - Discrete AHA", pos=(0.25, 0.95))

input("Press Enter to close QTriplot...")
q.close()