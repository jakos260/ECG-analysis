import numpy as np
import os, sys
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()
data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()  # add data path
root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))                         # add helpers

from qtripy import QTripy
from helpers.reader import read_geom_peacs_model
from helpers.helpers import load_multi_array

def get_test_pattern():
    # z
    # x y

    # matlab target
    # ver = [0 1 0; 2 0 0; 2 2 0; 1 1 2];
    # tri = [1 2 3; 1 4 2; 1 3 4; 2 4 3];

    vertices = np.array([
        [0.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 2.0, 0.0],    
        [1.0, 1.0, 2.0]
    ])
    triangles = np.array([
        [1,2,3],
        [1,4,2],
        [1,3,4],
        [2,4,3]
    ], dtype=np.int32)

    print("Vertices:")
    for i, v in enumerate(vertices, start=1):
        print(f"{i:2d}: {v}")

    return vertices, triangles


if __name__ == "__main__":

    # model = read_geom_peacs_model('./PPD1_ECGSIM_NEW','PPD1_ECGSIM_NEW')
    # vertices, triangles = model.GEOM.thorax.VER, model.GEOM.thorax.ITRI
    # vertices, triangles = model.GEOM.ventr.VER, model.GEOM.ventr.ITRI


    # refPos = model.GEOM.thorax.VER[38, :] - model.GEOM.thorax.NORMV[38, :] * 20
    # iatrial = 1006; atrialpos = model.ATRIA.geom.VER[iatrial-1,:] + model.ATRIA.geom.NORMV[iatrial-1,:] * 1
    # irvapex =   39; rvapexpos = model.GEOM.ventr.VER[irvapex-1,:] + model.GEOM.ventr.NORMV[irvapex-1,:] * 0.5
    # irvcoil =  119; rvcoilpos = model.GEOM.ventr.VER[irvcoil-1,:] + model.GEOM.ventr.NORMV[irvcoil-1,:] * 19
    # ilvprox =  962; lvproxpos = model.GEOM.ventr.VER[ilvprox-1,:] + model.GEOM.ventr.NORMV[ilvprox-1,:] * 0.5
    # ilvmid1 =  905; lvmid1pos = model.GEOM.ventr.VER[ilvmid1-1,:] + model.GEOM.ventr.NORMV[ilvmid1-1,:] * 0.5
    # ilvmid2 = 1207; lvmid2pos = model.GEOM.ventr.VER[ilvmid2-1,:] + model.GEOM.ventr.NORMV[ilvmid2-1,:] * 0.5
    # ilvdist =  636; lvdistpos = model.GEOM.ventr.VER[ilvdist-1,:] + model.GEOM.ventr.NORMV[ilvdist-1,:] * 0.5

    # pos_points=[rvapexpos,rvcoilpos,lvproxpos,lvmid1pos,lvmid2pos,lvdistpos,atrialpos]
    # for pos in pos_points:
    #     q.marker(pos,'red',5)

    # q.marker(refPos, 'black', 12)
    # q.surface(vertices, triangles, "demo")
    # q.cmd('color heart')
    # q.cmd('trans 0.6')
    # q.cmd('edge y')

    arrays = load_multi_array(os.path.join(data_path, "ECGsim_data", "PPD1_ECGSIM_NEW", "model", "ventricle.tri"))
    
    # arrays = load_multi_array(os.path.join(data_path, "ECGsim_data", "normal_young_male", "model", "ventricle.tri"))
    # rep = load_multi_array(os.path.join(data_path, "ECGsim_data", "normal_young_male", "ventricular_beats", "beat1", "user.rep"))[0]
    vertices, triangles = arrays[0][:,1:4], arrays[1][:,1:4]

    q = QTripy()
    q.begin()
    q.reset()
    q.surface(vertices, triangles)
    q.transparency(0.3)
    # q.values(rep, vmin=100)
    # q.gradient_bins(20)
    q.property_on_mouse_click('coor')
    q.text("ECGsim Ventricle", pos=(0.5, 0.95))
    # q.background_color("white")

    input("Press Enter to close QTriplot...")
    q.close()
