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
sample_data_path = os.path.join(data_path, "ECGsim_data", "normal_male")
# sample_data_path = os.path.join(data_path, "ECGsim_data", "normal_young_male")

heart_data = read_ecg_sim(sample_data_path)
ventricles_ver, ventriclec_tri = heart_data['VENTR']['geom']['VER'], heart_data['VENTR']['geom']['ITRI']
rcav_ver, rcav_tri = heart_data['GEOM']['rcav']['VER'], heart_data['GEOM']['rcav']['ITRI']
lcav_ver, lcav_tri = heart_data['GEOM']['lcav']['VER'], heart_data['GEOM']['lcav']['ITRI']

# if "normal_young_male" in sample_data_path:
#     shift = (-0.3702, -1.0695, 1.6489)
#     rcav_ver = rcav_ver + shift
#     lcav_ver = lcav_ver + shift

def intersect_rows(a, b):
    a_view = a.view([('', a.dtype)] * a.shape[1])
    b_view = b.view([('', b.dtype)] * b.shape[1])
    return np.intersect1d(a_view, b_view)

common_AB = len(intersect_rows(ventricles_ver, rcav_ver))
common_AC = len(intersect_rows(ventricles_ver, lcav_ver))

print("ventricles ∩ rcav:", common_AB)
print("ventricles ∩ lcav:", common_AC)

splitingMethod = SplitingMethod.SPLIT_BY_DISTANCE_TO_TYP
if splitingMethod == SplitingMethod.SPLIT_BY_DISTANCE:
    epicardium_ver, epicardium_tri, epicardium_ids, endocardium_ver, endocardium_tri, endocardium_ids = split(
        ventricles_ver,
        ventriclec_tri,
        rcav_ver,
        lcav_ver,
        method=splitingMethod,
        threshold=0.01
        )

    print(f"{max(endocardium_tri.flatten())=} {min(endocardium_tri.flatten())=} {endocardium_ver.shape=}")
    print(f"{max(epicardium_tri.flatten())=} {min(epicardium_tri.flatten())=} {epicardium_ver.shape=}")


    # save_tri_file(
    #     os.path.join(sample_data_path, "model/ventricle_endo.tri"),
    #     endocardium_ver,
    #     endocardium_tri,
    #     vertex_ids=endocardium_ids,
    # )
    # save_tri_file(
    #     os.path.join(sample_data_path, "model/ventricle_epi.tri"),
    #     epicardium_ver,
    #     epicardium_tri,
    #     vertex_ids=epicardium_ids,
    # )

    q = QTripy()
    q.begin()
    q.reset()

    q.markers(endocardium_ver, color='blue', r=1)
    # q.surface(endocardium_ver, endocardium_tri, color='blue')
    # q.edge('y')

    q.markers(epicardium_ver, color='red', r=1)
    # q.surface(epicardium_ver, epicardium_tri, color='red', opacity=0.5)
    # q.edge('y')

    q.markers(lcav_ver, color='black', r=1)
    q.markers(rcav_ver, color='black', r=1)
    # q.surface(lcav_ver, lcav_tri, color='black', opacity=0.7)
    # q.surface(rcav_ver, rcav_tri, color='black', opacity=0.7)
    # q.marker(np.mean(ventricles_ver, axis=0), 'white', 5)

    q.text("Epi and Endo split", pos=(0.25, 0.95))

    input("Press Enter to close QTriplot...")
    q.close()

if splitingMethod == SplitingMethod.SPLIT_BY_DISTANCE_TO_TYP:
    typ = split(
        ventricles_ver,
        ventriclec_tri,
        rcav_ver,
        lcav_ver,
        method=splitingMethod,
        threshold=0.01
        )
    
    with open(os.path.join(sample_data_path, "model/ventricle.typ"), 'w+') as f:
        for items in typ:
            f.write('%s\n' %items)

    print(f"done, file {os.path.join(sample_data_path, "model/ventricle.typ")} created")