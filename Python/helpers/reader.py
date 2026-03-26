"""
Python port of readGeomPeacsModel_prof.m (lightweight, numpy-based).
Provides read_geom_peacs_model_prof(dirname, subject, *args) -> MODEL
MODEL is a dynamic attribute container so you can do: MODEL.GEOM.thorax.VER
This port implements the main geometry/file readers and the helpers needed
for thorax triangle lookup used by AMA/lead mapping. It supports ASCII
matrix .lead/.aedl/.vedl files and simple .tri files. Binary .mat-like files
are only attempted as ASCII fallback; extend loadmat() if you need full
binary compatibility.
"""
from dataclasses import dataclass, field
import os
import numpy as np
import math
from typing import Any, Optional, Sequence
import io

# Small dynamic struct so nested fields can be created/accessed as attributes
class Struct:
    def __init__(self, d=None):
        if d:
            for k, v in d.items():
                setattr(self, k, v)
    def __repr__(self):
        keys = [k for k in self.__dict__.keys()]
        return f"Struct({keys})"
    def to_dict(self):
        return self.__dict__

def read_geom_peacs_model(dirname: str, subject: str, *args) -> Struct:
    """
    Usage:
      MODEL = read_geom_peacs_model('./PPD1_ECGSIM_NEW','PPD1_ECGSIM_NEW')
    Returns a Struct with nested attributes like MODEL.GEOM.thorax.VER, etc.
    """
    model = Struct()
    DATA = Struct()
    GEOM = Struct()
    DATA.GEOM = GEOM

    modeldir = os.path.join(dirname, subject) + '_'

    # helper to build file path and test existence
    def fpath(name):
        return modeldir + name

    # load atria if present
    if os.path.exists(fpath('atria.adj2d')):
        V, I = loadtri(fpath('atria.tri'))
        DATA.ATRIA = Struct()
        DATA.ATRIA.geom = Struct()
        DATA.ATRIA.geom.VER = V
        DATA.ATRIA.geom.ITRI = I
        NV, NT = trinormals(V, I)
        DATA.ATRIA.geom.NORMV = NV
        DATA.ATRIA.geom.NORMT = NT

        DATA.ATRIA.ADJsurf = loadmat(fpath('atria.adj2d'))
        DATA.ATRIA.ADJ3D = loadmat(fpath('atria.adj3d'))
        DATA.ATRIA.DISTsurf = loadmat(fpath('atria.dst2d'))
        DATA.ATRIA.DIST3D = loadmat(fpath('ATRIA.dst3d'))
        DATA.ATRIA.ADJANIS = loadmat(fpath('ATRIA.adjanis'))
        DATA.ATRIA.DISTANIS = loadmat(fpath('ATRIA.dstanis'))

        GEOM.atria = Struct()
        GEOM.atria.VER = V
        GEOM.atria.ITRI = I
        GEOM.atria.geomtyp = loadmat(fpath('atria.typ'))

    # ventricles (assume present)
    Vv, Iv = loadtri(fpath('ventricles.tri'))
    DATA.VENTR = Struct()
    DATA.VENTR.geom = Struct()
    DATA.VENTR.geom.VER = Vv
    DATA.VENTR.geom.ITRI = Iv
    NVv, NTv = trinormals(Vv, Iv)
    DATA.VENTR.geom.NORMV = NVv
    DATA.VENTR.geom.NORMT = NTv
    DATA.VENTR.ADJsurf = loadmat(fpath('ventricles.adj2d'))
    DATA.VENTR.ADJ3D = loadmat(fpath('ventricles.adj3d'))
    DATA.VENTR.DISTsurf = loadmat(fpath('ventricles.dst2d'))
    DATA.VENTR.DIST3D = loadmat(fpath('ventricles.dst3d'))
    DATA.VENTR.ADJANIS = loadmat(fpath('ventricles.adjanis'))
    DATA.VENTR.DISTANIS = loadmat(fpath('ventricles.dstanis'))
    DATA.VENTR.HEARTDIST = loadmat(fpath('ventricles.heartdist'))

    GEOM.ventr = Struct()
    GEOM.ventr.VER = Vv
    GEOM.ventr.ITRI = Iv
    NVg, NTg = trinormals(GEOM.ventr.VER, GEOM.ventr.ITRI)
    GEOM.ventr.NORMV = NVg
    GEOM.ventr.NORMT = NTg
    GEOM.ventr.geomtyp = loadmat(fpath('ventricles.typ'))
    GEOM.ventr.scar = loadmat(fpath('ventricles.scar'))
    GEOM.ventr.walls = loadmat(fpath('ventricles.walls'))
    GEOM.ventr.segments = loadmat(fpath('ventricles.segments'))

    # other geometry pieces (lcav, rcav, lungs, thorax, etc.)
    GEOM.lcav = Struct(); GEOM.lcav.VER, GEOM.lcav.ITRI = loadtri(fpath('lcav.tri'))
    GEOM.rcav = Struct(); GEOM.rcav.VER, GEOM.rcav.ITRI = loadtri(fpath('rcav.tri'))
    if os.path.exists(fpath('llung.tri')):
        GEOM.llung = Struct(); GEOM.llung.VER, GEOM.llung.ITRI = loadtri(fpath('llung.tri'))
    if os.path.exists(fpath('rlung.tri')):
        GEOM.rlung = Struct(); GEOM.rlung.VER, GEOM.rlung.ITRI = loadtri(fpath('rlung.tri'))

    GEOM.thorax = Struct()
    GEOM.thorax.VER, GEOM.thorax.ITRI = loadtri(fpath('thorax.tri'))
    GEOM.thorax.NORMV, GEOM.thorax.NORMT = trinormals(GEOM.thorax.VER, GEOM.thorax.ITRI)

    if os.path.exists(fpath('liver.tri')):
        GEOM.liver = Struct(); GEOM.liver.VER, GEOM.liver.ITRI = loadtri(fpath('liver.tri'))
    if os.path.exists(fpath('ribcage.tri')):
        GEOM.ribcage = Struct(); GEOM.ribcage.VER, GEOM.ribcage.ITRI = loadtri(fpath('ribcage.tri'))
    if os.path.exists(fpath('fatpad_1.tri')):
        GEOM.fatpad_1 = Struct(); GEOM.fatpad_1.VER, GEOM.fatpad_1.ITRI = loadtri(fpath('fatpad_1.tri'))
    if os.path.exists(fpath('fatpad_2.tri')):
        GEOM.fatpad_2 = Struct(); GEOM.fatpad_2.VER, GEOM.fatpad_2.ITRI = loadtri(fpath('fatpad_2.tri'))

    # handle lead files and wct detection (attempt to mimic MATLAB logic)
    wct = None
    # if third arg is path to leads
    if len(args) == 1 and isinstance(args[0], str):
        leads_dir = args[0]
    else:
        leads_dir = dirname

    # try find 'standard12lead' file in leads_dir
    lead_files = [f for f in os.listdir(leads_dir) if f.endswith('.lead')]
    std12 = None
    for fn in lead_files:
        if 'standard12lead' in fn:
            std12 = fn
            break
    if std12:
        ver = ensure_lead_xyz(loadmat(os.path.join(leads_dir, std12)))
        if ver.size:
            meanTh = np.mean(GEOM.thorax.VER, axis=0)
            wct_list = []
            for irow in range(min(3, ver.shape[0])):
                P = ver[irow, :]
                R = np.array([meanTh[0], meanTh[1], P[2]])
                TRIS = linetris(GEOM.thorax.VER, GEOM.thorax.ITRI, P, R)
                # filter by smallest abs distance
                if TRIS.size:
                    dists = np.abs(TRIS[:, -1])
                    TRIS = TRIS[dists <= dists.min()]
                    if TRIS.shape[0] >= 1:
                        tri_idx = int(TRIS[0, 0])
                        # take first vertex index of triangle (0-based)
                        wct_list.append(int(GEOM.thorax.ITRI[tri_idx, 0]))
            if wct_list:
                wct = wct_list
    # if no wct found and user passed an explicit wct (vector)
    if len(args) and not isinstance(args[0], str):
        wct = args[-1]

    # split AMA (attempt matching MAT logic; simplified)
    if wct is not None:
        DATA = splitAMA(modeldir, DATA, wct, True, GEOM)
        DATA = splitAMA(modeldir, DATA, wct, False, GEOM)

    # lead positions: set LEADPOS fields similar to MATLAB
    if std12:
        ver = ensure_lead_xyz(loadmat(os.path.join(leads_dir, std12)))
        if ver.size:
            name = 'standard12lead'
            iVer = np.zeros((len(ver,),), dtype=int)
            for i in range(len(ver)):
                d = np.linalg.norm(GEOM.thorax.VER - ver[i, :], axis=1)
                iVer[i] = int(np.argmin(d))
            if not hasattr(DATA.VENTR, 'LEADPOS'):
                DATA.VENTR.LEADPOS = Struct()
            DATA.VENTR.LEADPOS.__dict__[name] = ver
            DATA.VENTR.__dict__.setdefault('ILEADPOS', {})[name] = iVer

            if hasattr(DATA, 'ATRIA'):
                DATA.ATRIA.LEADPOS = Struct()
                DATA.ATRIA.LEADPOS.__dict__[name] = ver
                DATA.ATRIA.__dict__.setdefault('ILEADPOS', {})[name] = iVer
    else:
        # process all .lead files in dirname
        for fn in lead_files:
            ver = ensure_lead_xyz(loadmat(os.path.join(dirname, fn)))
            if ver.size:
                name = fn[:-5]
                name = name.replace(' ', '')
                name = name.replace(subject, '')
                name = 'A' + name
                name = name.replace('+', '_')
                AMA_A, AMA_V = getAMA(DATA, ver, GEOM)
                # store
                DATA.VENTR.__dict__[name] = AMA_V
                if not hasattr(DATA.VENTR, 'LEADPOS'):
                    DATA.VENTR.LEADPOS = Struct()
                DATA.VENTR.LEADPOS.__dict__[name] = ver
                if hasattr(DATA, 'ATRIA'):
                    DATA.ATRIA.__dict__[name] = AMA_A
                    if not hasattr(DATA.ATRIA, 'LEADPOS'):
                        DATA.ATRIA.LEADPOS = Struct()
                    DATA.ATRIA.LEADPOS.__dict__[name] = ver

    DATA.GEOM = GEOM
    # final small postprocessing to create ventr.endoVER and type approximations (kept simple)
    try:
        ventr_VER = GEOM.ventr.VER
        lcav_VER = GEOM.lcav.VER
        rcav_VER = GEOM.rcav.VER
        endoVER = np.zeros((len(ventr_VER),), dtype=int)
        for i in range(len(ventr_VER)):
            v = ventr_VER[i]
            if ((lcav_VER == v).all(axis=1)).any():
                endoVER[i] = 2
            if ((rcav_VER == v).all(axis=1)).any():
                endoVER[i] = 1
        GEOM.ventr.endoVER = endoVER
        GEOM.ventr.type = endoVER.copy()
    except Exception:
        pass

    model = DATA
    return model

# -------------------------
# Helper functions below
# -------------------------

def loadmat(name: str) -> np.ndarray:
    import numpy as np, os

    if not os.path.exists(name):
        raise FileNotFoundError(name)

    with open(name, "rb") as f:
        head8 = f.read(8)
        f.seek(0)

        # CASE A: specjalny nagłówek ';;mbfmat'
        if head8 == b';;mbfmat':
            f.read(1)  # skip 1 char
            _hs = int.from_bytes(f.read(4), 'little', signed=False)
            f.read(3)  # skip 3 chars
            rows = int.from_bytes(f.read(4), 'little', signed=False)
            cols = int.from_bytes(f.read(4), 'little', signed=False)
            need = rows * cols

            buf = f.read(need * 8)
            if len(buf) == need * 8:
                return np.frombuffer(buf, dtype='<f8', count=need).reshape((rows, cols))

            # fallback na float32
            # (pozycja po headerze: 8 + 1 + 4 + 3 + 8 = 24 bajty)
            f.seek(24)
            buf = f.read(need * 4)
            if len(buf) == need * 4:
                return np.frombuffer(buf, dtype='<f4', count=need).reshape((rows, cols)).astype(np.float64)
            raise ValueError("mbfmat: niepełne dane")

        # CASE B: surowy nagłówek (rows, cols) + binarne double/float
        rows_le = int.from_bytes(head8[:4], 'little', signed=False)
        cols_le = int.from_bytes(head8[4:8], 'little', signed=False)
        if rows_le > 0 and cols_le > 0:
            need = rows_le * cols_le
            rest = f.read()
            if len(rest) >= need * 8:
                return np.frombuffer(rest[:need*8], dtype='<f8', count=need).reshape((rows_le, cols_le))
            if len(rest) >= need * 4:
                return np.frombuffer(rest[:need*4], dtype='<f4', count=need).reshape((rows_le, cols_le)).astype(np.float64)

        # CASE C: ASCII — pierwsze 2 liczby to (rows cols), potem rows*cols wartości
        f.seek(0)
        txt = f.read().decode('utf-8', 'ignore').split()
        if len(txt) >= 2 and all(t.lstrip('-').isdigit() for t in txt[:2]):
            rows, cols = int(txt[0]), int(txt[1])
            need = rows * cols
            vals = np.array([float(x) for x in txt[2:2+need]], dtype=np.float64)
            if vals.size != need:
                raise ValueError(f"ASCII: oczekiwałem {need} liczb, mam {vals.size}")
            # Uwaga: w MATLAB jest fread([cols rows]) + M=M', co łącznie daje (rows, cols).
            return vals.reshape((rows, cols))

        # Ostateczny fallback (rzadko potrzebny)
        f.seek(0)
        arr = np.atleast_2d(np.loadtxt(f)).astype(np.float64)
        return arr

def loadtri(fn: str):
    """
    Robust line-oriented .tri reader:
      - first non-empty token: Npnt
      - next Npnt non-empty lines: "idx x y z"  (we take x,y,z)
      - next non-empty token: Ntri
      - next Ntri non-empty lines: "idx v1 v2 v3"  (we take v1..v3 and convert to 0-based)
    Returns (VER Nx3, ITRI Mx3) with 0-based indices. Tolerant to extra whitespace.
    """
    if not os.path.exists(fn):
        return np.array([]), np.array([], dtype=int)
    with open(fn, 'rt') as f:
        # collect non-empty stripped lines to keep parsing deterministic
        lines = [ln.strip() for ln in f if ln.strip() != '']
    if not lines:
        return np.array([]), np.array([], dtype=int)
    i = 0
    try:
        Npnt = int(lines[i].split()[0]); i += 1
    except Exception:
        raise RuntimeError(f"Unexpected .tri format (Npnt) in: {fn}")
    # read Npnt vertex lines
    if len(lines) < i + Npnt:
        raise RuntimeError(f"Unexpected .tri format: not enough vertex lines in {fn}")
    V = np.zeros((Npnt, 3), dtype=float)
    for k in range(Npnt):
        parts = lines[i + k].split()
        # typical: index x y z  -> take parts[1:4]; if not present try last 3 tokens
        if len(parts) >= 4:
            x, y, z = parts[1], parts[2], parts[3]
        elif len(parts) >= 3:
            # fallback (no leading index): take last three
            x, y, z = parts[-3], parts[-2], parts[-1]
        else:
            raise RuntimeError(f"Bad vertex line #{k+1} in {fn}: '{lines[i+k]}'")
        V[k, :] = [float(x), float(y), float(z)]
    i += Npnt
    # if no triangle block present -> return vertices only
    if i >= len(lines):
        return V, np.empty((0, 3), dtype=int)
    # read triangle count
    try:
        Ntri = int(lines[i].split()[0]); i += 1
    except Exception:
        # if cannot parse Ntri, return empty triangle list
        return V, np.empty((0, 3), dtype=int)
    if len(lines) < i + Ntri:
        # tolerate shorter triangle block by using available lines
        Ntri = max(0, len(lines) - i)
    ITRI = np.zeros((Ntri, 3), dtype=int)
    for k in range(Ntri):
        parts = lines[i + k].split()
        if len(parts) >= 4:
            # typical: idx v1 v2 v3
            v1, v2, v3 = parts[1], parts[2], parts[3]
        elif len(parts) >= 3:
            # fallback: last three tokens
            v1, v2, v3 = parts[-3], parts[-2], parts[-1]
        else:
            raise RuntimeError(f"Bad triangle line #{k+1} in {fn}: '{lines[i+k]}'")
        # convert to 0-based
        ITRI[k, :] = [int(v1) - 1, int(v2) - 1, int(v3) - 1]
    return V, ITRI

def trinormals(V: np.ndarray, ITRI: np.ndarray):
    """
    Match MATLAB trinormals.m implementation:
    NORMT(i,:) = cross(rp,rm)/2  with rp = v3-v1, rm = v2-v1
    NORMV at a vertex = sum( NORMT(tri)/3 ) over incident triangles
    finally normalize per-row.
    """
    if V is None or ITRI is None:
        return np.empty((0, 3)), np.empty((0, 3))

    tris = np.asarray(ITRI, dtype=int).copy()
    if tris.size == 0:
        return np.zeros((V.shape[0], 3)), np.empty((0, 3))

    if tris.ndim == 1:
        if tris.size % 3 == 0:
            tris = tris.reshape((-1, 3))
        else:
            return np.zeros((V.shape[0], 3)), np.empty((0, 3))

    ntri = tris.shape[0]
    nver = V.shape[0]

    # handle potential 1-based indices coming from other sources
    if tris.min() >= 1 and tris.max() <= nver:
        tris = tris - 1

    # gather triangle vertex coordinates
    v1 = V[tris[:, 0], :]
    v2 = V[tris[:, 1], :]
    v3 = V[tris[:, 2], :]

    # MATLAB: rp = v3-v1 ; rm = v2-v1 ; NORMT = cross(rp, rm)/2
    rp = v3 - v1
    rm = v2 - v1
    NORMT = np.cross(rp, rm) / 2.0

    # accumulate 1/3 of each triangle normal to its three vertices
    NORMV = np.zeros((nver, 3), dtype=float)
    thirds = NORMT / 3.0
    for i in range(ntri):
        i1, i2, i3 = tris[i, 0], tris[i, 1], tris[i, 2]
        n = thirds[i]
        NORMV[i1] += n
        NORMV[i2] += n
        NORMV[i3] += n

    # normalize rows, avoid divide-by-zero
    tnorms = np.linalg.norm(NORMT, axis=1)
    maskT = tnorms > 1e-12
    if np.any(maskT):
        NORMT[maskT] = NORMT[maskT] / tnorms[maskT, np.newaxis]

    vnorms = np.linalg.norm(NORMV, axis=1)
    maskV = vnorms > 1e-12
    if np.any(maskV):
        NORMV[maskV] = NORMV[maskV] / vnorms[maskV, np.newaxis]

    return NORMV, NORMT

def linetris(ver: np.ndarray, itri: np.ndarray, P: Sequence[float], R: Sequence[float]):
    """
    Intersect ray R + t*(P-R) with each triangle; return array of rows:
    [tri_idx, 0, u, v, t] where u,v are parameters in X = v1 + u*(v2-v1) + v*(v3-v1)
    and t is the ray parameter. Only triangles with valid barycentric coords are returned.
    Sorted by increasing abs(t).
    """
    P = np.asarray(P, dtype=float)
    R = np.asarray(R, dtype=float)
    d = P - R
    tris = itri.astype(int)
    rows = []
    eps = 1e-9
    for i in range(tris.shape[0]):
        v1 = ver[tris[i, 0], :]
        v2 = ver[tris[i, 1], :]
        v3 = ver[tris[i, 2], :]
        n = np.cross(v2 - v1, v3 - v1)
        denom = np.dot(n, d)
        if abs(denom) < 1e-12:
            continue
        t = np.dot(n, (v1 - R)) / denom
        X = R + t * d
        # solve [v2-v1, v3-v1] * [u;v] = X-v1
        M = np.column_stack((v2 - v1, v3 - v1))
        try:
            uv = np.linalg.lstsq(M, (X - v1), rcond=None)[0]
        except Exception:
            continue
        u, v = uv[0], uv[1]
        if u >= -1e-6 and v >= -1e-6 and (u + v) <= 1 + 1e-6:
            rows.append([i, 0.0, u, v, t])
    if len(rows) == 0:
        return np.empty((0, 5))
    arr = np.array(rows)
    # sort by abs(t)
    order = np.argsort(np.abs(arr[:, -1]))
    return arr[order, :]

def doWCT(Ain: np.ndarray, Awct: np.ndarray):
    if Ain is None or Ain.size == 0:
        return Ain
    return Ain - np.ones((Ain.shape[0], 1)) @ Awct.reshape((1, -1))

def calcAwct(AthorsoIn: np.ndarray, wct):
    if AthorsoIn is None or AthorsoIn.size == 0:
        return np.zeros((1,))
    wct = list(wct)
    return np.mean(AthorsoIn[wct, :], axis=0)

def ensure_2d(mat):
    """
    Ensure numeric input is a 2-D numpy array.
    - empty -> empty (0,0)
    - 0-D -> empty (0,0)
    - 1-D -> reshape to (-1,1)
    - 2-D -> return as float copy
    - higher dims -> flatten then reshape to (-1,1)
    """
    if mat is None:
        return np.empty((0, 0), dtype=float)
    a = np.asarray(mat)
    if a.size == 0:
        return np.empty((0, 0), dtype=float)
    if a.ndim == 0:
        return np.empty((0, 0), dtype=float)
    if a.ndim == 1:
        return a.reshape((-1, 1)).astype(float)
    if a.ndim >= 2:
        return a.astype(float)
    return a.astype(float)

def splitAMA(modeldir: str, DATA: Struct, wct, usetriag: bool, GEOM: Struct):
    """
    Simplified splitAMA: read *all.* or thorax.* aedl/vedl and slice into blocks
    following the MATLAB order. This version closely follows the MATLAB code:
    sequentially consume rows for THORAX, RCAV, LCAV, RLUNG, LLUNG, LIVER,
    RIBCAGE, FATPAD_1, FATPAD_2 and finally VENTRICLES (or tVENTRICLES).
    """
    extA = 'aedl'; extV = 'vedl'
    if usetriag:
        extA = 't' + extA; extV = 't' + extV

    # load and scale (40x) if present
    AA = None
    AV = None
    if os.path.exists(modeldir + 'all.' + extA):
        AA = loadmat(modeldir + 'all.' + extA)
    elif os.path.exists(modeldir + 'thorax.' + extA):
        AA = loadmat(modeldir + 'thorax.' + extA)
    if AA is not None:
        AA = 40.0 * AA

    if os.path.exists(modeldir + 'all.' + extV):
        AV = loadmat(modeldir + 'all.' + extV)
    elif os.path.exists(modeldir + 'thorax.' + extV):
        AV = loadmat(modeldir + 'thorax.' + extV)
    if AV is not None:
        AV = 40.0 * AV

    # ensure 2D arrays for consistent slicing
    AA = ensure_2d(AA)
    AV = ensure_2d(AV)

    # helper to safely set attribute on nested structs
    def _set(dest, name, val):
        try:
            setattr(dest, name, val)
        except Exception:
            dest.__dict__[name] = val

    # process ventricular blocks (AV)
    if isinstance(AV, np.ndarray) and AV.size:
        n = len(GEOM.thorax.VER)
        try:
            Awct = calcAwct(AV[:n, :], wct)
            AV = doWCT(AV, Awct)
        except Exception:
            pass

        prefix = 't' if usetriag else ''
        _set(DATA.VENTR, prefix + 'THORAX', AV[0:n, :])
        # advance index
        idx = n

        # sequentially try to take each block if present
        def _take_blk_geom(attrname, destname):
            nonlocal idx
            if hasattr(GEOM, attrname):
                length = len(getattr(GEOM, attrname).VER)
                if AV.shape[0] > idx:
                    end = min(idx + length, AV.shape[0])
                    _set(DATA.VENTR, prefix + destname, AV[idx:end, :])
                idx += length

        # follow MATLAB order
        _take_blk_geom('rcav', 'RCAV')
        _take_blk_geom('lcav', 'LCAV')
        if hasattr(GEOM, 'rlung'):
            _take_blk_geom('rlung', 'RLUNG')
        if hasattr(GEOM, 'llung'):
            _take_blk_geom('llung', 'LLUNG')
            if hasattr(GEOM, 'liver'):
                _take_blk_geom('liver', 'LIVER')
            if hasattr(GEOM, 'ribcage'):
                _take_blk_geom('ribcage', 'RIBCAGE')
            if hasattr(GEOM, 'fatpad_1'):
                _take_blk_geom('fatpad_1', 'FATPAD_1')
            if hasattr(GEOM, 'fatpad_2'):
                _take_blk_geom('fatpad_2', 'FATPAD_2')

        # final ventricles block if remaining rows equal ventr.VER length
        vlen = len(GEOM.ventr.VER)
        rem = AV.shape[0] - idx
        if rem == vlen:
            # take last vlen rows (MATLAB uses end-length+1:end)
            block = np.asarray(AV[-vlen:, :])
            # ensure 2-D
            if block.ndim == 1:
                block = block.reshape((-1, 1))
            # if shape doesn't have vlen rows but vlen columns -> transpose
            if block.shape[0] != vlen and block.shape[1] == vlen:
                block = block.T
            _set(DATA.VENTR, prefix + 'VENTRICLES', block)
        else:
            # fallback: if there are enough rows at idx to slice vlen, prefer that
            if AV.shape[0] >= idx + vlen:
                block = np.asarray(AV[idx: idx + vlen, :])
                if block.ndim == 1:
                    block = block.reshape((-1, 1))
                if block.shape[0] != vlen and block.shape[1] == vlen:
                    block = block.T
                _set(DATA.VENTR, prefix + 'VENTRICLES', block)

    # process atrial blocks (AA) similarly
    if isinstance(AA, np.ndarray) and AA.size:
        n = len(GEOM.thorax.VER)
        try:
            Awct = calcAwct(AA[:n, :], wct)
            AA = doWCT(AA, Awct)
        except Exception:
            pass

        prefix = 't' if usetriag else ''
        _set(DATA.ATRIA, prefix + 'THORAX', AA[0:n, :])
        idx = n

        def _take_blk_geom_A(attrname, destname):
            nonlocal idx
            if hasattr(GEOM, attrname):
                length = len(getattr(GEOM, attrname).VER)
                if AA.shape[0] > idx:
                    end = min(idx + length, AA.shape[0])
                    _set(DATA.ATRIA, prefix + destname, AA[idx:end, :])
                idx += length

        _take_blk_geom_A('lcav', 'LCAV')
        _take_blk_geom_A('rcav', 'RCAV')
        if hasattr(GEOM, 'rlung'):
            _take_blk_geom_A('rlung', 'RLUNG')
        if hasattr(GEOM, 'llung'):
            _take_blk_geom_A('llung', 'LLUNG')
        if hasattr(GEOM, 'ribcage'):
            _take_blk_geom_A('ribcage', 'RIBCAGE')
        if hasattr(GEOM, 'fatpad_1'):
            _take_blk_geom_A('fatpad_1', 'FATPAD_1')
        if hasattr(GEOM, 'fatpad_2'):
            _take_blk_geom_A('fatpad_2', 'FATPAD_2')

        # final atria block if present
        if hasattr(GEOM, 'atria'):
            alen = len(GEOM.atria.VER)
            rem = AA.shape[0] - idx
            if rem == alen:
                block = np.asarray(AA[-alen:, :])
                if block.ndim == 1:
                    block = block.reshape((-1, 1))
                if block.shape[0] != alen and block.shape[1] == alen:
                    block = block.T
                _set(DATA.ATRIA, prefix + 'ATRIA', block)
            else:
                if AA.shape[0] >= idx + alen:
                    block = np.asarray(AA[idx: idx + alen, :])
                    if block.ndim == 1:
                        block = block.reshape((-1, 1))
                    if block.shape[0] != alen and block.shape[1] == alen:
                        block = block.T
                    _set(DATA.ATRIA, prefix + 'ATRIA', block)

    return DATA

def getAMA(DATA: Struct, ver: np.ndarray, GEOM: Struct):
    """
    Compute AMA_A, AMA_V: interpolated thorax potentials at positions 'ver'.
    This implementation uses barycentric interpolation on thorax triangles.
    Returns (AMA_A, AMA_V) where AMA_A may be empty if DATA.ATRIA not present.
    """
    meanTh = np.mean(GEOM.thorax.VER, axis=0)
    AMA_A = None
    if hasattr(DATA, 'ATRIA') and hasattr(DATA.ATRIA, 'THORAX'):
        AMA_A = np.zeros((ver.shape[0], DATA.ATRIA.THORAX.shape[1]))
    AMA_V = np.zeros((ver.shape[0], DATA.VENTR.THORAX.shape[1]))

    for i in range(ver.shape[0]):
        P = ver[i, :]
        R = np.array([meanTh[0], meanTh[1], P[2]])
        TRIS = linetris(GEOM.thorax.VER, GEOM.thorax.ITRI, P, R)
        if TRIS.size == 0:
            # fallback: nearest vertex
            d = np.linalg.norm(GEOM.thorax.VER - P, axis=1)
            vidx = int(np.argmin(d))
            if AMA_A is not None:
                AMA_A[i, :] = DATA.ATRIA.THORAX[vidx, :]
            AMA_V[i, :] = DATA.VENTR.THORAX[vidx, :]
            continue
        # pick first intersection
        itri_idx = int(TRIS[0, 0])
        u = float(TRIS[0, 2]); v = float(TRIS[0, 3])
        itri = GEOM.thorax.ITRI[itri_idx, :]
        i1, i2, i3 = itri[0], itri[1], itri[2]
        w1 = 1.0
        # MATLAB formula used: value = val(i1) + (val(i2)-val(i1))*u + (val(i3)-val(i1))*v
        if AMA_A is not None:
            AMA_A[i, :] = DATA.ATRIA.THORAX[i1, :] + (DATA.ATRIA.THORAX[i2, :] - DATA.ATRIA.THORAX[i1, :]) * u + (DATA.ATRIA.THORAX[i3, :] - DATA.ATRIA.THORAX[i1, :]) * v
        AMA_V[i, :] = DATA.VENTR.THORAX[i1, :] + (DATA.VENTR.THORAX[i2, :] - DATA.VENTR.THORAX[i1, :]) * u + (DATA.VENTR.THORAX[i3, :] - DATA.VENTR.THORAX[i1, :]) * v

    return AMA_A, AMA_V

def ensure_lead_xyz(ver):
    """
    Normalize a loaded lead array into shape (N,3) containing X,Y,Z coordinates.
    Accepts 1D or 2D numeric arrays that may represent Nx3 or Nx4 (index+xyz)
    or a flattened sequence of coordinates.
    Returns an empty (0,3) array if data cannot be interpreted.
    """
    ver = np.asarray(ver)
    if ver.size == 0:
        return np.empty((0, 3), dtype=float)
    # 1D: try reshape to (N,4) -> take cols 1:4, else (N,3)
    if ver.ndim == 1:
        if ver.size % 4 == 0:
            ver = ver.reshape((-1, 4))[:, 1:4]
        elif ver.size % 3 == 0:
            ver = ver.reshape((-1, 3))
        else:
            # try to recover by truncating to multiple of 3
            n3 = (ver.size // 3) * 3
            if n3 == 0:
                return np.empty((0, 3), dtype=float)
            ver = ver[:n3].reshape((-1, 3))
    else:
        # 2D: common cases: Nx4 (index + xyz) or Nx3
        if ver.shape[1] >= 4:
            ver = ver[:, 1:4]
        elif ver.shape[1] == 3:
            ver = ver.copy()
        elif ver.shape[0] == 3 and ver.shape[1] != 3:
            # maybe transposed; try transpose
            vt = ver.T
            if vt.shape[1] == 3:
                ver = vt.copy()
            else:
                flat = ver.flatten()
                if flat.size % 3 == 0:
                    ver = flat.reshape((-1, 3))
                else:
                    return np.empty((0, 3), dtype=float)
        else:
            # fallback: try flatten & reshape to Nx3
            flat = ver.flatten()
            if flat.size % 3 == 0:
                ver = flat.reshape((-1, 3))
            else:
                return np.empty((0, 3), dtype=float)
    return ver.astype(float)
