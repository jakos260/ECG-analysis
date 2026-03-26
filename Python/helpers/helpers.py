import numpy as np

def gets(T, dep, rep, p, mode):
    """
    Python port of gets.m
    T: (nn, nt) array (time grid rows = vertices, columns = time samples)
    dep, rep: (nn,) arrays of times (can be column or row vectors)
    p: sequence of parameters (at least length 4; p[4] optional)
    mode: if ==1 use simplified formula
    Returns S: (nn, nt) array
    """
    T = np.asarray(T, dtype=float)
    dep = np.asarray(dep, dtype=float).reshape(-1)
    rep = np.asarray(rep, dtype=float).reshape(-1)
    p = np.asarray(p, dtype=float).reshape(-1)

    if T.shape[0] != dep.shape[0]:
        # try transpose if user provided T as (nt,nn)
        if T.shape[1] == dep.shape[0]:
            T = T.T
        else:
            raise ValueError("T and dep size mismatch")

    nn, nt = T.shape

    # broadcasted time differences
    TDEP = T - dep[:, None]

    p1 = float(p[0]) if p.size > 0 else 1.0
    p2 = float(p[1]) if p.size > 1 else 0.0
    p3 = float(p[2]) if p.size > 2 else 1.0
    p4 = float(p[3]) if p.size > 3 else 1.0
    p5 = float(p[4]) if p.size > 4 else 0.0

    if mode != 1:
        TREP = T - rep[:, None] - p5
        # Y = (p2 + 1/(1+exp(p3*TREP))) / (1+exp(p4*TREP))
        with np.errstate(over='ignore', invalid='ignore'):
            denom1 = 1.0 + np.exp(p3 * TREP)
            term1 = 1.0 / denom1
            denom2 = 1.0 + np.exp(p4 * TREP)
            Y = (p2 + term1) / denom2

        with np.errstate(over='ignore', invalid='ignore'):
            S = (1.0 - Y) / (1.0 + np.exp(-p1 * TDEP))

        # normalize each row by its max (re-establish unit upstroke)
        row_max = np.nanmax(S, axis=1)
        # avoid division by zero
        row_max[row_max == 0] = 1.0
        S = S / row_max[:, None]
    else:
        with np.errstate(over='ignore', invalid='ignore'):
            S = 1.0 / (1.0 + np.exp(-p1 * TDEP))

    # replace any non-finite with nan
    S = np.asarray(S, dtype=np.float64)
    S[~np.isfinite(S)] = np.nan
    return S


def getTris_prof(dep, rep, pp, maxt, VER, ITRI):
    """
    Python port of getTris_prof.m
    dep, rep: arrays of vertex times (length = nver)
    pp: parameter vector where pp[2], pp[3] correspond to MATLAB pp(3), pp(4)
    maxt: integer, number of time steps (MATLAB used 1..maxt)
    VER: (nver,3) vertex coords
    ITRI: (ntri,3) triangle indices (1-based or 0-based). Function accepts either.
    Returns S: (ntri, maxt) array of triangle activation fractions / TMP
    """
    VER = np.asarray(VER, dtype=float)
    ITRI = np.asarray(ITRI, dtype=int).copy()
    dep = np.asarray(dep, dtype=float).reshape(-1)
    rep = np.asarray(rep, dtype=float).reshape(-1)
    pp = np.asarray(pp, dtype=float).reshape(-1)

    # handle 1-based indices
    if ITRI.size > 0 and ITRI.min() == 1:
        ITRI = ITRI - 1

    ntri = ITRI.shape[0]
    S = np.zeros((ntri, int(maxt)), dtype=float)

    # compute triangle areas
    v0 = VER[ITRI[:, 0], :]
    v1 = VER[ITRI[:, 1], :]
    v2 = VER[ITRI[:, 2], :]
    crossp = np.cross(v1 - v0, v2 - v0)
    triArea = 0.5 * np.linalg.norm(crossp, axis=1)
    # avoid division by zero later
    triArea_safe = triArea.copy()
    triArea_safe[triArea_safe == 0] = 1.0

    # parameters for rep processing (MATLAB pp(3), pp(4) are pp[2], pp[3] here)
    pp3 = float(pp[2]) if pp.size > 2 else 1.0
    pp4 = float(pp[3]) if pp.size > 3 else 1.0

    # iterate faces (keeps logic identical to MATLAB)
    for iFace in range(ntri):
        idx = ITRI[iFace]
        vv0 = VER[idx[0], :]
        vv1 = VER[idx[1], :]
        vv2 = VER[idx[2], :]

        tdep0 = dep[idx[0]]
        tdep1 = dep[idx[1]]
        tdep2 = dep[idx[2]]

        area = triArea[iFace]
        if area == 0.0:
            # degenerate triangle -> keep zeros
            continue

        for t in range(1, int(maxt) + 1):  # MATLAB 1..maxt
            # follow the same boolean checks as MATLAB code
            if (t >= tdep0) or (t > tdep1) or (t >= tdep2):
                # determine which vertex leads (smallest dep)
                if (tdep0 <= tdep1) and (tdep0 <= tdep2):
                    # vertex 0 leads
                    # compute fractions safely
                    f01 = 0.0
                    denom = (tdep1 - tdep0)
                    if denom == 0:
                        f01 = 1.0 if t >= tdep1 else 0.0
                    else:
                        f01 = np.clip((t - tdep0) / denom, 0.0, 1.0)

                    f02 = 0.0
                    denom = (tdep2 - tdep0)
                    if denom == 0:
                        f02 = 1.0 if t >= tdep2 else 0.0
                    else:
                        f02 = np.clip((t - tdep0) / denom, 0.0, 1.0)

                    v01 = (vv1 - vv0) * f01
                    v02 = (vv2 - vv0) * f02
                    triangleFraction = np.linalg.norm(np.cross(v01, v02)) / 2.0
                    S[iFace, t - 1] = triangleFraction / area

                elif tdep1 <= tdep2:
                    # vertex 1 leads
                    denom = (tdep0 - tdep1)
                    if denom == 0:
                        f10 = 1.0 if t >= tdep0 else 0.0
                    else:
                        f10 = np.clip((t - tdep1) / denom, 0.0, 1.0)

                    denom = (tdep2 - tdep1)
                    if denom == 0:
                        f12 = 1.0 if t >= tdep2 else 0.0
                    else:
                        f12 = np.clip((t - tdep1) / denom, 0.0, 1.0)

                    v10 = (vv0 - vv1) * f10
                    v12 = (vv2 - vv1) * f12
                    triangleFraction = np.linalg.norm(np.cross(v10, v12)) / 2.0
                    S[iFace, t - 1] = triangleFraction / area

                else:
                    # vertex 2 leads
                    denom = (tdep0 - tdep2)
                    if denom == 0:
                        f20 = 1.0 if t >= tdep0 else 0.0
                    else:
                        f20 = np.clip((t - tdep2) / denom, 0.0, 1.0)

                    denom = (tdep1 - tdep2)
                    if denom == 0:
                        f21 = 1.0 if t >= tdep1 else 0.0
                    else:
                        f21 = np.clip((t - tdep2) / denom, 0.0, 1.0)

                    v20 = (vv0 - vv2) * f20
                    v21 = (vv1 - vv2) * f21
                    triangleFraction = np.linalg.norm(np.cross(v20, v21)) / 2.0
                    S[iFace, t - 1] = triangleFraction / area

                if S[iFace, t - 1] > 0.999:
                    trepstart = t
                    S[iFace, t - 1] = 1.0
                    # if rep provided, compute trailing shape for t >= trepstart
                    if rep is not None and rep.size > 0:
                        trep0 = rep[idx[0]]
                        trep1 = rep[idx[1]]
                        trep2 = rep[idx[2]]
                        for tt in range(trepstart, int(maxt) + 1):
                            # compute s as sum of two logistic terms per vertex (MATLAB expression)
                            term0 = (1.0 / (1.0 + np.exp(pp3 * (tt - trep0)))) * (1.0 / (1.0 + np.exp(pp4 * (tt - trep0))))
                            term1 = (1.0 / (1.0 + np.exp(pp3 * (tt - trep1)))) * (1.0 / (1.0 + np.exp(pp4 * (tt - trep1))))
                            term2 = (1.0 / (1.0 + np.exp(pp3 * (tt - trep2)))) * (1.0 / (1.0 + np.exp(pp4 * (tt - trep2))))
                            s = term0 + term1 + term2
                            S[iFace, tt - 1] = 1.0 - (s / 3.0)
                    break

    return S

def zeromean(a, axis=None, nan_safe=False):
    """
    Return array with mean along `axis` removed.
    - a: array-like
    - axis: None (global mean), 0 (column mean), 1 (row mean), etc.
    - nan_safe: if True use nanmean to ignore NaNs
    """
    a = np.asarray(a, dtype=float)
    if nan_safe:
        m = np.nanmean(a, axis=axis, keepdims=True)
    else:
        m = np.mean(a, axis=axis, keepdims=True)
    return a - m


def load_multi_array(path):
    """Load a text file that contains one or more arrays concatenated.

    Format expected:
      - A block begins with a line containing an integer N (number of rows).
      - The next N non-empty, non-comment lines are rows of that array (space separated values).
      - After those N rows another integer may appear, starting the next block, and so on until EOF.

    Returns a list of numpy arrays. For each block, the function detects the
    maximum number of columns among the block's rows and returns a 2D numpy
    array of shape (N, max_cols). Rows with fewer columns are padded with NaN.
    Comment lines starting with `#` or `%` and empty lines are ignored.
    """
    arrays = []
    with open(path, 'r', encoding='utf-8') as fh:
        raw_lines = [ln.rstrip('\n') for ln in fh]

    # filter out empty lines and comments
    lines = []
    for ln in raw_lines:
        s = ln.strip()
        if not s:
            continue
        if s[0] in ('#', '%'):
            continue
        lines.append(s)

    i = 0
    while i < len(lines):
        # parse block row count
        parts = lines[i].split()
        if not parts:
            i += 1
            continue
        try:
            # allow counts like '10' or '10.0' by casting through float
            count = int(float(parts[0]))
        except Exception as e:
            raise ValueError(f"Expected integer row count at line {i+1}: '{lines[i]}'") from e
        i += 1

        # read the following `count` rows
        block_rows = []
        for r in range(count):
            if i >= len(lines):
                raise ValueError(f"Unexpected end of file while reading {count} rows (needed row {r+1}).")
            row_tokens = lines[i].split()
            try:
                nums = [float(x) for x in row_tokens]
            except Exception:
                # if conversion fails produce empty row
                nums = []
            block_rows.append(nums)
            i += 1

        # determine max columns and build array, padding with NaN where needed
        maxcols = max((len(r) for r in block_rows), default=0)
        if maxcols == 0:
            arr = np.empty((len(block_rows), 0), dtype=float)
        else:
            arr = np.full((len(block_rows), maxcols), np.nan, dtype=float)
            for ridx, row in enumerate(block_rows):
                if row:
                    arr[ridx, :len(row)] = row

        arrays.append(arr)

    return arrays