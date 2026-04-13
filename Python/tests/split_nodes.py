import os, sys
import itertools
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
from qtripy.readECGsim import read_ecg_sim
from helpers.save_tri_file import save_tri_file
from helpers.split_endo_and_epi_nodes import SplitingMethod, split


load_dotenv()
data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
sample_data_path = os.path.join(data_path, "ECGsim_data", "normal_young_male")

heart_data = read_ecg_sim(sample_data_path)
ventricles_ver, ventriclec_tri = heart_data['VENTR']['geom']['VER'], heart_data['VENTR']['geom']['ITRI']
rcav_ver, rcav_tri = heart_data['GEOM']['rcav']['VER'], heart_data['GEOM']['rcav']['ITRI']
lcav_ver, lcav_tri = heart_data['GEOM']['lcav']['VER'], heart_data['GEOM']['lcav']['ITRI']

# Old approach is preserved for reference.
# epicardium_ver, epicardium_tri, endocardium_ver, endocardium_tri = split_vertices_by_distance(
#     ventricles_ver, ventriclec_tri, rcav_ver, lcav_ver, threshold=2.0)

# New method: initial distance split plus cleanup of small disconnected triangle groups.
epicardium_ver, epicardium_tri, epicardium_ids, endocardium_ver, endocardium_tri, endocardium_ids = split(
    ventricles_ver,
    ventriclec_tri,
    rcav_ver,
    lcav_ver,
    method=SplitingMethod.SPLIT_BY_DISTANCE_WITH_COMPONENT_CLEANUP,
    threshold=2.0,
    min_component_triangles=40,
    min_component_share=0.005,
    rcav_tri=rcav_tri,
    lcav_tri=lcav_tri,
)

print(f"{max(endocardium_tri.flatten())=} {min(endocardium_tri.flatten())=} {endocardium_ver.shape=}")
print(f"{max(epicardium_tri.flatten())=} {min(epicardium_tri.flatten())=} {epicardium_ver.shape=}")

input_ids = np.arange(1, len(ventricles_ver) + 1, dtype=int)
all_split_ids = np.unique(np.concatenate([epicardium_ids, endocardium_ids]))
overlap_ids = np.intersect1d(epicardium_ids, endocardium_ids)
missing_ids = np.setdiff1d(input_ids, all_split_ids)
extra_ids = np.setdiff1d(all_split_ids, input_ids)
print(f"Input vertices: {len(input_ids)}")
print(f"Epicardium vertices: {len(epicardium_ids)} | Endocardium vertices: {len(endocardium_ids)}")
print(f"Unique split IDs: {len(all_split_ids)} | Overlap IDs: {len(overlap_ids)}")
print(f"Missing input IDs: {len(missing_ids)} | Extra split IDs: {len(extra_ids)}")
if len(missing_ids) == 0 and len(extra_ids) == 0:
    print("Split coverage test: PASS")
    if len(overlap_ids) > 0:
        print(f"Shared boundary vertex IDs in both outputs: {len(overlap_ids)}")
else:
    print("Split coverage test: FAIL")
    if len(missing_ids) > 0:
        print(f"Missing vertex IDs from split: {missing_ids[:20]}{'...' if len(missing_ids) > 20 else ''}")
    if len(extra_ids) > 0:
        print(f"Extra vertex IDs in split: {extra_ids[:20]}{'...' if len(extra_ids) > 20 else ''}")
    if len(overlap_ids) > 0:
        print(f"Overlapping epicardium/endocardium IDs: {len(overlap_ids)} IDs shared between both outputs")

save_tri_file(
    os.path.join(sample_data_path, "model/ventricle_endo.tri"),
    endocardium_ver,
    endocardium_tri,
    vertex_ids=endocardium_ids,
)
save_tri_file(
    os.path.join(sample_data_path, "model/ventricle_epi.tri"),
    epicardium_ver,
    epicardium_tri,
    vertex_ids=epicardium_ids,
)

q = QTripy()
q.begin()
q.reset()

q.markers(endocardium_ver, color='blue', r=1)
q.surface(endocardium_ver, endocardium_tri, color='blue')
q.edge('y')

q.markers(epicardium_ver, color='red', r=1)
q.surface(epicardium_ver, epicardium_tri, color='red', opacity=0.5)
q.edge('y')

# q.markers(lcav_ver, color='black', r=1)
# q.markers(rcav_ver, color='black', r=1)
# q.marker(np.mean(ventricles_ver, axis=0), 'white', 5)

q.text("Epi and Endo split", pos=(0.25, 0.95))

input("Press Enter to close QTriplot...")
q.close()