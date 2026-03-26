import os
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from dotenv import load_dotenv

# optional - for signal smoothing puposes
# from scipy.ndimage import gaussian_filter1d


def read_tri(file_path):
    with open(file_path, 'r') as f:
        n_points = int(f.readline())
        points = np.loadtxt(f, max_rows=n_points).reshape((n_points, 4))[:, 1:]

        try:
            n_triangles = int(f.readline())
            triangles = np.loadtxt(f, max_rows=n_triangles).reshape((n_triangles, 4))[:, 1:] - 1
        except Exception:
            triangles = np.array([])
    return points, triangles


def loadmat_py(filepath):
    if not os.path.exists(filepath):
        return np.array([]), ""

    try:
        with open(filepath, 'r') as f:
            header = f.readline().split()
            if len(header) == 2:
                N = list(map(int, header))
                M = np.loadtxt(f, max_rows=N[0])
                if M.shape != (N[0], N[1]):
                    M = M.reshape((N[0], N[1]))
                return M.T, ""  
    except Exception:
        pass 

   
    with open(filepath, 'rb') as f:
        magic = f.read(8)
        f.seek(0)  

        if magic == b';;mbfmat':
            f.read(1)  
            f.read(4)  
            f.read(1)
            f.read(1)
            f.read(1)
            N = np.fromfile(f, dtype=np.int32, count=2)
            M = np.fromfile(f, dtype=np.float64, count=N[0] * N[1])
            M = M.reshape((N[1], N[0])).T
        else:
            f.seek(0)
            N = np.fromfile(f, dtype=np.int32, count=2)
            M = np.fromfile(f, dtype=np.float32, count=N[0] * N[1])
            M = M.reshape((N[1], N[0])).T

        # Read extra info at the end
        extra = f.read(1000)
        try:
            extra = extra.decode('ascii').strip()
        except UnicodeDecodeError:
            extra = ''

    return M, extra


def load_lead12(file_path):
    M, extraresults = loadmat_py(file_path)
    return np.r_[M[0:3], M[9:12], M[3:9]]



def calc_awct(AthorsoIn, wct):
    Awct = np.mean(AthorsoIn[wct, :], axis=0)
    return Awct


def do_wct(Ain, Awct):
    if Ain.size > 0:
        return Ain - np.ones((Ain.shape[0], 1)) * Awct
    else:
        return Ain



def getS(T, dep, rep, p, mode=4):
    """
    Python version of MATLAB getS function
    T    : time matrix (nn x nt)
    dep  : depolarization times (nn x 1)
    rep  : repolarization times (nn x 1)
    p    : parameter vector [p1..p5]
    mode : behavior switch (default 4)
    """
    nn, nt = T.shape

    dep = dep.reshape((nn, 1))
    rep = rep.reshape((nn, 1))

    # Shift times
    TDEP = T - dep
    TREP = T - rep - p[4]

    Y = (p[1] + 1.0 / (1.0 + np.exp(p[2] * TREP))) / (1.0 + np.exp(p[3] * TREP))

    S = (1.0 - Y) / (1.0 + np.exp(-p[0] * TDEP))

    S = S / np.max(S, axis=1, keepdims=True)

    return S


def read_ecg_sim(dirname):
    DATA = {}

    modeldir = os.path.join(dirname, 'model')
    ecgdir = os.path.join(dirname, 'ecgs')
    abeatsdir = os.path.join(dirname, 'atrial_beats', 'beat1')
    vbeatsdir = os.path.join(dirname, 'ventricular_beats', 'beat1')

    def fexists(path):
        return os.path.isfile(path)

    atria_tri = os.path.join(modeldir, 'atria.tri')
    if fexists(atria_tri):
        DATA['ATRIA'] = {}
        VER, ITRI = read_tri(atria_tri)
        geom = {'VER': VER, 'ITRI': ITRI}
        atria = DATA['ATRIA']
        atria['geom'] = geom
        atria['ADJsurf'] = loadmat_py(os.path.join(modeldir, 'atria.adj2d'))[0]
        atria['ADJ3D'] = loadmat_py(os.path.join(modeldir, 'atria.adj3d'))[0]
        atria['DISTsurf'] = loadmat_py(os.path.join(modeldir, 'atria.dst2d'))[0]
        atria['DIST3D'] = loadmat_py(os.path.join(modeldir, 'ATRIA.dst3d'))[0]
        atria['ADJANIS'] = loadmat_py(os.path.join(modeldir, 'ATRIA.adjanis'))[0]
        atria['DISTANIS'] = loadmat_py(os.path.join(modeldir, 'ATRIA.dstanis'))[0]
        atria['THORAX'] = loadmat_py(os.path.join(modeldir, 'atria2Thorax.mat'))[0]
        atria['ATRIA'] = loadmat_py(os.path.join(modeldir, 'atria2atria.mat'))[0]
        atria['VENTRICLES'] = loadmat_py(os.path.join(modeldir, 'atria2Ventricles.mat'))[0]
        atria['RLUNG'] = loadmat_py(os.path.join(modeldir, 'atria2RLung.mat'))[0]
        atria['LLUNG'] = loadmat_py(os.path.join(modeldir, 'atria2LLung.mat'))[0]
        atria['LCAV'] = loadmat_py(os.path.join(modeldir, 'atria2LCavity.mat'))[0]
        atria['RCAV'] = loadmat_py(os.path.join(modeldir, 'atria2RCavity.mat'))[0]

        try:
            atria['lead12'] = loadmat_py(os.path.join(modeldir, 'atria2standard12lead.mat'))[0]
        except:
            atria['lead12'] = loadmat_py(os.path.join(modeldir, 'atria2standard_12.mat'))[0]
            atria['VENTR212lead'] = atria['lead12'][np.r_[0:3, 9:12, 3:9], :]

    DATA['VENTR'] = {}
    ventr = DATA['VENTR']
    ventr['geom'] = {}
    ventr['geom']['VER'], ventr['geom']['ITRI'] = read_tri(os.path.join(modeldir, 'ventricle.tri'))
    ventr['ADJsurf'] = loadmat_py(os.path.join(modeldir, 'ventricle.adj2d'))[0]
    ventr['ADJ3D'] = loadmat_py(os.path.join(modeldir, 'ventricle.adj3d'))[0]
    ventr['DISTsurf'] = loadmat_py(os.path.join(modeldir, 'ventricle.dst2d'))[0]
    ventr['DIST3D'] = loadmat_py(os.path.join(modeldir, 'ventricle.dst3d'))[0]
    ventr['ADJANIS'] = loadmat_py(os.path.join(modeldir, 'ventricle.adjanis'))[0]
    ventr['DISTANIS'] = loadmat_py(os.path.join(modeldir, 'ventricle.dstanis'))[0]
    ventr['THORAX'] = loadmat_py(os.path.join(modeldir, 'ventricles2Thorax.mat'))[0]
    ventr['ATRIA'] = loadmat_py(os.path.join(modeldir, 'ventricles2atria.mat'))[0]
    ventr['VENTRICLES'] = loadmat_py(os.path.join(modeldir, 'ventricles2Ventricles.mat'))[0]
    ventr['RLUNG'] = loadmat_py(os.path.join(modeldir, 'ventricles2RLung.mat'))[0]
    ventr['LLUNG'] = loadmat_py(os.path.join(modeldir, 'ventricles2LLung.mat'))[0]
    ventr['LCAV'] = loadmat_py(os.path.join(modeldir, 'ventricles2LCavity.mat'))[0]
    ventr['RCAV'] = loadmat_py(os.path.join(modeldir, 'ventricles2RCavity.mat'))[0]

    try:
        ventr['VENTR212lead'] = loadmat_py(os.path.join(modeldir, 'ventricles2standard12lead.mat'))[0]
    except:
        ventr['lead12'] = loadmat_py(os.path.join(modeldir, 'ventricles2standard_12.mat'))[0]
        ventr['lead12'] = ventr['lead12'][np.r_[0:3, 9:12, 3:9], :]

    # Geometry loading
    DATA['GEOM'] = {}
    geom = DATA['GEOM']
    if fexists(atria_tri):
        geom['atria'] = {}
        geom['atria']['VER'], geom['atria']['ITRI'] = read_tri(atria_tri)

    for part in ['ventricle', 'lcav', 'rcav', 'llung', 'rlung', 'thorax']:
        geom[part] = {}
        geom[part]['VER'], geom[part]['ITRI'] = read_tri(os.path.join(modeldir, f'{part}.tri'))

    # Alias ventr to ventricle
    geom['ventr'] = geom['ventricle']

    # Compute endoVER
    geom['ventr']['endoVER'] = np.zeros(len(geom['ventr']['VER']))
    for i, v in enumerate(geom['ventr']['VER']):
        if any(np.all(v == geom['lcav']['VER'], axis=1)):
            geom['ventr']['endoVER'][i] = 2
        elif any(np.all(v == geom['rcav']['VER'], axis=1)):
            geom['ventr']['endoVER'][i] = 1

    if 'atria' in geom:
        geom['atria']['endoVER'] = np.zeros(len(geom['atria']['VER']))
        for i, v in enumerate(geom['atria']['VER']):
            if any(np.all(v == geom['lcav']['VER'], axis=1)):
                geom['atria']['endoVER'][i] = 2
            elif any(np.all(v == geom['rcav']['VER'], axis=1)):
                geom['atria']['endoVER'][i] = 1

    # Load ECG files
    DATA['ECG'] = {}
    for file in os.listdir(ecgdir):
        if file.endswith('.refECG'):
            key = file[:-7].replace('(', '').replace(')', '').replace(' ', '_')
            elec_key = key + 'elec'

            if 'standard_12' in key:
                value = loadmat_py(os.path.join(ecgdir, file[:-7]) + '.elec')[0]
                DATA["ECG"][elec_key] = value

                if value.shape[0] == 9 and abs(value[-1, 3] - value[-2, 3]) > 200:
                    DATA["ECG"][key] = load_lead12(os.path.join(ecgdir, file))
                    value = np.vstack([value[6:], value[:6]])
                    DATA["ECG"][elec_key] = value

            else:
                DATA["ECG"][key] = loadmat_py(os.path.join(ecgdir, file))[0]
                DATA["ECG"][elec_key] = loadmat_py(os.path.join(ecgdir, file[:-7]) + '.elec')[0]
            
            if DATA["ECG"][elec_key] is None or len(DATA["ECG"][elec_key]) == 0:
                DATA["ECG"][elec_key] = DATA["GEOM"]["thorax"]["VER"]
            else:
                 DATA["ECG"][elec_key] = DATA["ECG"][elec_key][:, 1:4]


    # Load atrial and ventricular beats
    DATA['ABEAT'] = {}
    if os.path.exists(abeatsdir) and os.path.isdir(abeatsdir):
        for file in os.listdir(abeatsdir):
            if file.startswith('user.'):
                key = file[5:]
                DATA['ABEAT'][key] = loadmat_py(os.path.join(abeatsdir, file))[0]

    DATA['VBEAT'] = {}
    if os.path.exists(vbeatsdir) and os.path.isdir(vbeatsdir):
        for file in os.listdir(vbeatsdir):
            if file.startswith('user.'):
                key = file[5:]
                DATA['VBEAT'][key] = loadmat_py(os.path.join(vbeatsdir, file))[0]

    # WCT lead detection
    ECGelec = DATA['ECG'].get('standard12leadelec') or DATA['ECG'].get('standard_12elec')
    index = np.zeros(len(ECGelec), dtype=int)
    for i, lead in enumerate(ECGelec):
        dist = np.linalg.norm(geom['thorax']['VER'] - lead, axis=1)
        min_idx = np.argmin(dist)
        if dist[min_idx] < 20:
            index[i] = min_idx

    wct = index[:3]
    geom['wct'] = wct
    geom['stand12leadsIndex'] = index

    # WCT subtraction
    if 'ATRIA' in DATA and DATA['ATRIA']['THORAX'].shape[0] > 1:
        Awct = calc_awct(DATA['ATRIA']['THORAX'], wct)
        for field in ['THORAX', 'ATRIA', 'VENTRICLES', 'RLUNG', 'LLUNG', 'RCAV', 'LCAV']:
            DATA['ATRIA'][field] = do_wct(DATA['ATRIA'][field], Awct)

    if DATA['VENTR']['THORAX'].shape[0] > 1:
        Awct = calc_awct(DATA['VENTR']['THORAX'], wct)
        for field in ['THORAX', 'ATRIA', 'VENTRICLES', 'RLUNG', 'LLUNG', 'RCAV', 'LCAV']:
            DATA['VENTR'][field] = do_wct(DATA['VENTR'][field], Awct)

    return DATA


def simulate_ecg_signal(DATA, node_index=19, maxt=550, foci=[100], fociact=[0]):
    """
    Generate ECG signal at given thorax node using depolarization and repolarization.
    node_index : index of thorax node to extract ECG signal from (0-based)
    maxt       : max simulation time in ms
    foci       : list of foci node indices
    fociact    : activation offsets for each focus
    """
    # Load heart-to-thorax transfer matrix
    HD = DATA['VENTR']['THORAX']  # nn x nt

    # Step 1: depolarization times (per heart node)
    dep = np.zeros(HD.shape[1])  # nn = number of heart nodes (columns)
    for f, t0 in zip(foci, fociact):
        dep += HD[f - 1, :] + t0
    dep /= len(foci)
    dep = dep.reshape(-1, 1)  # shape (nn, 1)

    # Step 2: repolarization times
    rep = 250 + np.mean(dep) - dep * 0.6

    # Step 3: time vector
    T = np.tile(np.arange(0, maxt + 1), (len(dep), 1))

    # Step 4: TMP simulation
    p = np.array([2, 0, -0.025, -0.02, 0])
    S = getS(T, dep, rep, p, mode=4)
    # S = gaussian_filter1d(S, sigma=2, axis=1)  # optional smoothing

    # Step 5: ECG lead extraction
    ECG = HD[node_index - 1, :] @ S  
    ECG /= np.max(np.abs(ECG))       

    return ECG, np.arange(0, maxt + 1)


if __name__ == "__main__":
    load_dotenv()
    data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
    
    my_data = read_ecg_sim(os.path.join(data_path, "ECGsim_data", "normal_male"))
    ECG_signal, time_vector = simulate_ecg_signal(
        my_data,
        node_index=19,
        maxt=550,
        foci=[100],
        fociact=[0]
    )

    plt.figure(figsize=(10, 4))
    plt.plot(time_vector, ECG_signal, 'r')
    plt.title("Simulated ECG signal at thorax node 19")
    plt.xlabel("Time (ms)")
    plt.ylabel("Amplitude (a.u.)")
    plt.grid(True)
    plt.show()