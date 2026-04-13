import os
import socket
import struct
import subprocess
import time
import numpy as np
from pathlib import Path
from dotenv import load_dotenv

class QTripyCmd:
    def __send(self, cmd):
        self.cmd(cmd)

    def transparency(self, alpha):
        self.__send(f"trans {alpha}")

    def gradient_step(self, step):
        self.__send(f"step {step}")

    def reset(self):
        self.__send("reset")

    def edge(self, color):
        self.__send(f"edge {color}")

    def property_on_mouse_click(self, property_name):
        '''coor, tri, ver, fun, lm, cross'''
        self.__send(f"mouse {property_name}")

    def background_color(self, color):
        self.__send(f"bgdcolor {color}")

    def color_range(self, vmin, vmax):
        self.minmax_values = (vmin, vmax)
        self.__send(f"funscale {vmin} {vmax}")

    def text(self, text, pos=(0,0)):
        self.__send(f'textcolor black')
        if any(x>1 or x<0 for x in pos):
            raise ValueError("pos must be in (0,1) range for both x and y")
        cmd = f'text {pos[0]} {pos[1]} {text}'
        self.__send(cmd)


class QTripy(QTripyCmd):
    def __init__(self, path = None):
        self.qtriplot_path = path
        if self.qtriplot_path == None:
            load_dotenv()
            self.qtriplot_path = Path(os.getenv("ENV_QTRIPLOT_PATH")).resolve()
        if not os.path.exists(self.qtriplot_path):
            raise FileNotFoundError(f"The system cannot find the {self.qtriplot_path}")

        self.socket = None
        self._proc = None
        self._marker_counter = 0
        self._objects = []
        self.panels = (1,1) # horizontal x vertical
        self.active_panel = (1,1)
        self.minmax_values = None
        self.min_time_before_next_cmd = 0.0  # debouncing: minimum time between surface/values commands in seconds
        self.last_cmd_time = 0.0  # timestamp of the last surface/values command

    def __start_qtriplot(self, port):
        for host in ["127.0.0.1", "::1", "localhost"]:
            try:
                s = socket.create_connection((host, port), timeout=1)
                print(f"Connected to QTriplot at {host}:{port}")
                return s
            except Exception:
                continue

        self._proc = subprocess.Popen([self.qtriplot_path, "-p", str(port)])

        for _ in range(50):
            time.sleep(0.2)
            for host in ["127.0.0.1", "::1", "localhost"]:
                try:
                    s = socket.create_connection((host, port), timeout=1)
                    print(f"Connected to QTriplot at {host}:{port}")
                    return s
                except Exception:
                    continue
        raise RuntimeError("Cannot connect to QTriplot")

    def begin(self, port=1041):
        self.socket = self.__start_qtriplot(port)
        self.__reset_panels()

    def cmd(self, cmd, enable_print=False):
        msg = cmd.encode("ascii")
        self.socket.sendall(struct.pack(">i", len(msg)))
        self.socket.sendall(msg)
        if enable_print:
            print("Command:", cmd)

    def __should_debounce_cmd(self):
        """
        Check if a surface/values command should be debounced.
        Returns True if the command should be skipped due to debouncing.
        """
        if self.min_time_before_next_cmd <= 0:
            return False
        
        current_time = time.time()
        time_since_last = current_time - self.last_cmd_time
        
        if time_since_last < self.min_time_before_next_cmd:
            return True
        
        self.last_cmd_time = current_time
        return False

    def enable_debounce(self, seconds):
        """Enable debouncing with the given minimum interval in seconds.

        Passing a non-positive value disables debouncing (equivalent to
        calling `disable_debounce`).
        """
        try:
            sec = float(seconds)
        except Exception:
            sec = 0.0
        if sec <= 0.0:
            self.min_time_before_next_cmd = 0.0
            return
        self.min_time_before_next_cmd = sec

    def disable_debounce(self):
        """Disable debouncing (allow all commands through)."""
        self.min_time_before_next_cmd = 0.0

    def _place_marker(self, pos, color, r, name=None, set_active=True, set_color=True):
        VER, ITRI = self.__make_sphere(nref=3)
        P = np.asarray(pos, dtype=float).reshape(1, 3)
        verts = VER * float(r) + P

        self.surface(verts.astype(np.float64), ITRI.astype(np.int32), name=name)

        if set_color:
            self.cmd(f"color {color}")

        # optionally update active panel once, not for every point in batch mode
        if set_active:
            self.__set_active_panel()

    def marker(self, pos, color, r, name=None):
        self._place_marker(pos, color, r, name=name)

    def markers(self, positions, color, r, names=None, set_active_after=True, combine=True):
        """Place multiple markers in one batch.

        positions: iterable of (x,y,z) positions
        color: color name or RGB string for all markers
        r: radius for marker spheres
        names: optional iterable of names, one per marker
        set_active_after: if True, activates panel once at end (best for bulk updates)
        combine: if True, send as one combined surface command for all markers
        cut_points_outside_enable: if True, drop markers outside self.max_distance
        """
        positions = list(positions)
        if names is not None:
            names = list(names)
            if len(names) != len(positions):
                raise ValueError("names must have same length as positions")

        if combine and positions:
            vertices_list = []
            triangles_list = []
            for i, pos in enumerate(positions):
                name = names[i] if names is not None else None
                VER, ITRI = self.__make_sphere(nref=3)
                P = np.asarray(pos, dtype=float).reshape(1, 3)
                vertices_list.append(VER * float(r) + P)
                triangles_list.append(ITRI)

            # one combined surface command
            self.surface(vertices_list, triangles_list, name=name)
            self.cmd(f"color {color}")

        else:
            for i, pos in enumerate(positions):
                name = names[i] if names is not None else None
                self._place_marker(pos, color, r, name=name, set_active=False)

        if set_active_after:
            self.__set_active_panel()

    def create_cylinder(self, length, radius=3, nseg=32, nlen=1, axis=(0, 0, 1), cap_ends=True, name="cylinder"):
        VER, ITRI, CNTR_IDX = self.__make_cylinder(length, radius, nseg, nlen, axis, cap_ends)
        self.surface(VER.astype(np.float64), ITRI.astype(np.int32), name=name)
        return VER, ITRI, CNTR_IDX

    def _combine_meshes(self, vertices_list, triangles_list):
        """Concatenate multiple mesh parts into one mesh."""
        if not vertices_list or not triangles_list:
            return np.empty((0, 3), dtype=float), np.empty((0, 3), dtype=int)

        v_acc = []
        t_acc = []
        offset = 0
        for v, t in zip(vertices_list, triangles_list):
            v = np.asarray(v, dtype=float)
            t = np.asarray(t, dtype=int)
            v_acc.append(v)
            t_acc.append(t + offset)
            offset += v.shape[0]

        all_verts = np.vstack(v_acc) if v_acc else np.empty((0, 3), dtype=float)
        all_tris = np.vstack(t_acc) if t_acc else np.empty((0, 3), dtype=int)
        return all_verts, all_tris

    def surface(self, vertices, triangles, name=None, enable_print=False, color=None, opacity=None):
        # Accept either a single mesh or lists of meshes
        if isinstance(vertices, (list, tuple)) and isinstance(triangles, (list, tuple)):
            vertices, triangles = self._combine_meshes(vertices, triangles)

        # Prepare the numeric mesh before computing distances
        vertices = np.asarray(vertices, dtype=float)
        triangles = np.asarray(triangles, dtype=int)

        # Check debouncing
        if self.__should_debounce_cmd():
            return
        
        nver = vertices.shape[0]
        ntri = triangles.shape[0] if triangles.size > 0 else 0
        object_name = f"obj_{self._marker_counter}" if name is None else f"{name}_{self._marker_counter}"
        self._marker_counter += 1
        # record created object for later deletion
        try:
            self._objects.append(object_name)
        except Exception:
            # defensive: if list can't be appended for some reason, ignore
            pass
        header = f"|tri {object_name}\0"
        while len(header) % 8 != 0:
            header += "\0"
        ver_data = vertices.T.astype(">f8").tobytes()
        tri_data = triangles.astype(">i4").T.tobytes() if ntri > 0 else b""
        tot = len(header) + 8 + 8 * nver * 3 + 4 * ntri * 3
        self.socket.sendall(struct.pack(">i", tot))
        self.socket.sendall(header.encode("ascii"))
        self.socket.sendall(struct.pack(">ii", nver, ntri))
        self.socket.sendall(ver_data)
        self.socket.sendall(tri_data)
        if enable_print:
            print(f"Surface: {nver} verts, {ntri} tris")
        if opacity is not None:
            self.transparency(opacity)
        if color:
            self.cmd(f"color {color}")


    def values(self, fun, name='', vmin=None, vmax=None):
        """
        Send function values to qtriplot to color the surface.
        Args:
            fun: numpy array of values (N,) or (N,M)
            name: optional name of the geometry to apply values to
            auto_scale: if True (default) compute vmin/vmax from data and send range command
            vmin, vmax: optional explicit range (overrides auto computed values)
        """
        # Check debouncing
        if self.__should_debounce_cmd():
            return
        
        fun = np.asarray(fun, dtype=float)
        if fun.ndim == 1:
            nr, nc = fun.shape[0], 1
        else:
            nr, nc = fun.shape

        # determine range
        dmin = float(np.nanmin(fun)) if fun.size else 0.0
        dmax = float(np.nanmax(fun)) if fun.size else 1.0

        if vmin:
            dmin = vmin
        if vmax:
            dmax = vmax
            
        if self.minmax_values is None or vmin or vmax:
            self.color_range(dmin, dmax)

        # if vmin is None or vmax is None:
        #     if auto_scale:
        #         # compute data range ignoring NaNs
        #         dmin = float(np.nanmin(fun)) if fun.size else 0.0
        #         dmax = float(np.nanmax(fun)) if fun.size else 1.0
        #         # if constant array, expand small epsilon to avoid zero range
        #         if dmax <= dmin:
        #             eps = max(1.0, abs(dmin)) * 1e-6
        #             dmin -= eps
        #             dmax += eps
        #         vmin = dmin if vmin is None else vmin
        #         vmax = dmax if vmax is None else vmax
        #     else:
        #         # fallback defaults (preserve previous behaviour)
        #         vmin = 0.0 if vmin is None else vmin
        #         vmax = float(np.nanmax(fun)) if vmax is None else vmax

        # # send a range command so qtriplot can use correct color mapping
        # try:
        #     self.cmd(f"fun range {name} {vmin} {vmax}")
        # except Exception:
        #     # non-fatal: continue even if qtriplot doesn't support the command
        #     pass

        # Format command string with padding to 8-byte boundary
        command = f"|fun {name}\0"
        while len(command) % 8 != 0:
            command += "\0"

        # Calculate total message size: command + 2 int32s + double array
        tot = len(command) + 8 + 8 * nr * nc

        # Send header
        self.socket.sendall(struct.pack(">i", tot))
        self.socket.sendall(command.encode("ascii"))
        
        # Send dimensions
        self.socket.sendall(struct.pack(">ii", nr, nc))
        
        # Send values in big-endian double format
        self.socket.sendall(fun.astype(">f8").tobytes())

        # ensure active panel is set after sending values
        self.__set_active_panel()


    def gradient_bins(self, bins):
        if bins < 1:
            bins = 10
        step = (self.minmax_values[1] - self.minmax_values[0]) / bins
        return super().gradient_step(step)

    
    def set_active_panel(self, horizontal=1, vertical=1):
        """
        Set the active panel in QTriplot.
        Args:
            horizontal: horizontal panel index (1-based)
            vertical: vertical panel index (1-based)
        """
        self.active_panel = (min(horizontal, self.panels[0]), min(vertical, self.panels[1]))

    def __set_active_panel(self):
        if self.socket:
            self.cmd(f'panel {self.active_panel[0]} {self.active_panel[1]}')


    def set_panels_number(self, horizontal=1, vertical=1):
        """
        Set the number of horizontal and vertical panels in QTriplot.
        Args:
            horizontal: number of horizontal panels (default 1)
            vertical: number of vertical panels (default 1)
        """
        self.panels = (horizontal, vertical)
        if self.socket:
            self.cmd(f'horizontal {horizontal}')
            self.cmd(f'vertical {vertical}')

    def __reset_panels(self):
        """
        Reset the panel layout to a single panel (1x1).
        """
        self.set_panels_number(horizontal=1, vertical=1)



    def __make_sphere(self, nref=3):
        """
        Build an icosphere (unit radius) and return (VER, ITRI)
        VER: (N,3) float numpy array
        ITRI: (M,3) int numpy array (0-based indices)
        nref: number of midpoint refinements (>=0)
        """
        # create base icosahedron
        t = (1.0 + np.sqrt(5.0)) / 2.0
        verts = np.array([
            [-1,  t,  0],
            [ 1,  t,  0],
            [-1, -t,  0],
            [ 1, -t,  0],
            [ 0, -1,  t],
            [ 0,  1,  t],
            [ 0, -1, -t],
            [ 0,  1, -t],
            [ t,  0, -1],
            [ t,  0,  1],
            [-t,  0, -1],
            [-t,  0,  1],
        ], dtype=float)

        # normalize to unit sphere
        verts /= np.linalg.norm(verts, axis=1)[:, None]

        faces = np.array([
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ], dtype=int)

        # refine
        for _ in range(int(max(0, round(nref)))):
            midpoint_cache = {}
            new_faces = []
            verts_list = verts.tolist()

            def get_midpoint(i1, i2):
                key = tuple(sorted((int(i1), int(i2))))
                if key in midpoint_cache:
                    return midpoint_cache[key]
                v1 = np.array(verts_list[key[0]])
                v2 = np.array(verts_list[key[1]])
                m = (v1 + v2) / 2.0
                m /= np.linalg.norm(m)
                verts_list.append(m.tolist())
                idx = len(verts_list) - 1
                midpoint_cache[key] = idx
                return idx

            for tri in faces:
                a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
                ba = get_midpoint(b, a)
                ac = get_midpoint(a, c)
                cb = get_midpoint(c, b)
                new_faces.extend([
                    [b, ba, cb],
                    [ba, a, ac],
                    [cb, ac, c],
                    [ba, ac, cb],
                ])
            verts = np.array(verts_list, dtype=float)
            faces = np.array(new_faces, dtype=int)

        # ensure unit length
        verts /= np.linalg.norm(verts, axis=1)[:, None]

        return verts.astype(float), faces.astype(int)

    def __make_cylinder(self, length, radius=3, nseg=32, nlen=1, axis=(0, 0, 1), cap_ends=True):
        """
        Build a cylinder and return (VER, ITRI).
        length: total cylinder length
        radius: cylinder radius
        nseg: number of segments around the circumference (>=3)
        nlen: number of segments along the length (>=1)
        axis: orientation vector (length-3); cylinder is centered at origin along this axis
        cap_ends: whether to add flat end caps
        Returns:
          VER: (N,3) float numpy array
          ITRI: (M,3) int numpy array (0-based indices)
        """

        axis = np.asarray(axis, dtype=float).reshape(3)
        norm = np.linalg.norm(axis)
        if norm == 0:
            raise ValueError("axis must be non-zero")
        axis_u = axis / norm

        # build orthonormal basis (u,v,axis_u)
        if abs(axis_u[0]) < 0.9:
            tmp = np.array([1.0, 0.0, 0.0])
        else:
            tmp = np.array([0.0, 1.0, 0.0])
        u = np.cross(axis_u, tmp)
        u /= np.linalg.norm(u)
        v = np.cross(axis_u, u)
        v /= np.linalg.norm(v)

        # ring positions along axis (centered). we create nlen segments -> nlen+1 rings
        rings = int(max(1, int(nlen))) + 1
        zs = np.linspace(-length / 2.0, length / 2.0, rings)

        verts_list = []
        seg = int(max(3, int(nseg)))
        for z in zs:
            center = axis_u * float(z)
            for k in range(seg):
                theta = 2.0 * np.pi * k / float(seg)
                pos = center + radius * (np.cos(theta) * u + np.sin(theta) * v)
                verts_list.append(pos.tolist())

        verts = np.asarray(verts_list, dtype=float)
        faces = []

        # side faces (quads -> two tris per quad)
        for i in range(rings - 1):
            base = i * seg
            base_next = (i + 1) * seg
            for j in range(seg):
                jn = (j + 1) % seg
                v0 = base + j
                v1 = base + jn
                v2 = base_next + j
                v3 = base_next + jn
                # keep consistent winding so normals point outward
                faces.append([v0, v2, v1])
                faces.append([v1, v2, v3])

        # caps
        if cap_ends:
            # bottom center vertex (at zs[0])
            bottom_center_idx = verts.shape[0]
            verts = np.vstack([verts, (axis_u * zs[0]).reshape(1, 3)])
            # top center vertex (at zs[-1])
            top_center_idx = verts.shape[0]
            verts = np.vstack([verts, (axis_u * zs[-1]).reshape(1, 3)])

            # bottom cap: ring 0. Want normals pointing outward (i.e. -axis_u for bottom if axis points +)
            base = 0
            for j in range(seg):
                jn = (j + 1) % seg
                # triangle orientation: (center, next, current) -> outward normal
                faces.append([base + jn, bottom_center_idx, base + j])

            # top cap: last ring. Want normals pointing outward (+axis_u)
            base = (rings - 1) * seg
            for j in range(seg):
                jn = (j + 1) % seg
                # triangle orientation: (center, current, next) -> outward normal
                faces.append([base + j, top_center_idx, base + jn])

        verts = verts.astype(float)
        faces = np.asarray(faces, dtype=int)

        return verts, faces, (bottom_center_idx, top_center_idx)

    def close(self):
        self.clear()        
        if self._proc:
            self._proc.terminate()
            self._proc.wait()
            print("QTriplot process terminated.")
        if self.socket:
            self.socket.close()
            self.socket = None
            print("Connection closed.")

    def axis(self, x_len=10.0, y_len=25.0, z_len=50.0, radius=1.0, segments=32, offset=(0.0, -100.0, -50.0)):
        """Create and display three positive-axis cylinders (X, Y, Z).

        The cylinders are aligned with the positive X, Y and Z axes, all
        starting at the origin and extending in the positive direction.

        Args:
            x_len: length of the cylinder along the X axis.
            y_len: length of the cylinder along the Y axis.
            z_len: length of the cylinder along the Z axis.
            radius: radius of all cylinders.
            segments: number of radial subdivisions for each cylinder.
            offset: 3D offset added to all cylinder vertices.
        """
        offset = np.asarray(offset, dtype=float)
        if offset.shape != (3,):
            raise ValueError('offset must be a 3-element vector')

        # Create cylinders along each axis
        # Cylinders start at origin and extend in positive direction
        axes_config = [
            (x_len, (1.0, 0.0, 0.0), (x_len / 2.0, 0.0, 0.0), 'x-axis'),
            (y_len, (0.0, 1.0, 0.0), (0.0, y_len / 2.0, 0.0), 'y-axis'),
            (z_len, (0.0, 0.0, 1.0), (0.0, 0.0, z_len / 2.0), 'z-axis'),
        ]

        for length, axis_vec, pos_offset, name in axes_config:
            verts, tris, _ = self.__make_cylinder(
                length=float(length),
                radius=float(radius),
                nseg=int(segments),
                nlen=1,
                axis=axis_vec,
                cap_ends=True
            )
            # Shift cylinder from centered to positive direction, then apply user offset
            verts = verts + np.asarray(pos_offset, dtype=float) + offset
            self.surface(verts.astype(np.float64), tris.astype(np.int32), name=name)
            self.cmd(f"color black")

    def clear(self, enable_print=False):
        """
        Delete all objects previously sent to qtriplot and reset internal list/counter.
        """
        if not self.socket:
            # nothing to do if not connected
            self._objects = []
            self._marker_counter = 0
            return

        # try deleting each recorded object (best-effort)
        for obj in list(self._objects):
            try:
                # primary delete command
                self.cmd(f"delete {obj}", enable_print=enable_print)
            except Exception:
                try:
                    # fallback variant
                    self.cmd(f"del {obj}", enable_print=enable_print)
                except Exception:
                    # ignore errors - best-effort cleanup
                    pass

        # reset bookkeeping
        self._objects = []
        self._marker_counter = 0





