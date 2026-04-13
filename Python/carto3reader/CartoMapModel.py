import os
import re
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET

class CartoPoint:
    def __init__(self, point_id, xml_path, base_dir):
        self.id = point_id
        self.xml_path = xml_path
        self.base_dir = base_dir
        
        self.bipolar = None
        self.unipolar = None
        self.ecg_filename = None
        self.pos_filename = None
        
        self._position_df = None
        self._ecg_df = None
        
        self.__parse_point_xml()

    def __parse_point_xml(self):
        tree = ET.parse(self.xml_path)
        root = tree.getroot()
        
        v_node = root.find('Voltages')
        if v_node is not None:
            self.bipolar = float(v_node.get('Bipolar', 0))
            self.unipolar = float(v_node.get('Unipolar', 0))
            
        ecg_node = root.find('ECG')
        if ecg_node is not None:
            self.ecg_filename = ecg_node.get('FileName')
            
        for conn in root.findall('.//Connector'):
            for attr_name, attr_val in conn.attrib.items():
                if 'OnAnnotation' in attr_val:
                    self.pos_filename = attr_val
                    break

    def __read_carto_tsv(self, path, **kwargs):
        header = 0
        try:
            with open(path, 'r', encoding='utf-8', errors='replace') as f:
                first = f.readline()
                second = f.readline()
        except FileNotFoundError:
            raise

        if first and '\t' not in first and second and '\t' in second:
            header = 1

        return pd.read_csv(path, sep='\t', header=header, index_col=False, **kwargs)

    @property
    def position(self):
        if self._position_df is None and self.pos_filename:
            path = os.path.join(self.base_dir, self.pos_filename)
            if os.path.exists(path):
                df = self.__read_carto_tsv(path, usecols=['X', 'Y', 'Z'], nrows=1)
                if not df.empty:
                    self._position_df = df.iloc[0].values
        return self._position_df

    def ecg(self):
        if self._ecg_df is None and self.ecg_filename:
            path = os.path.join(self.base_dir, self.ecg_filename)
            if os.path.exists(path):
                self._ecg_df = self.__read_carto_tsv(path)
        return self._ecg_df


class CartoMapModel:
    def __init__(self):
        self.points = []

        self.mesh_vertices = None
        self.mesh_triangles = None

        self.mesh_colors = None
        self.color_names = []
        
        self.max_distance = None  # threshold for filtering vertices

    def _load_master(self):
        tree = ET.parse(self.master_xml)
        for p_entry in tree.getroot().findall('Point'):
            p_id = p_entry.get('ID')
            p_xml = os.path.join(self.base_dir, p_entry.get('File_Name'))

            if os.path.exists(p_xml):
                self.points.append(CartoPoint(p_id, p_xml, self.base_dir))

    def load_mesh(self, mesh_filename=None):
        """Loads heart anatomy and associated color data from .mesh file"""
        if mesh_filename is None:
            mesh_path = self.mesh_path
        else:
            mesh_path = os.path.join(self.base_dir, mesh_filename)

        if not os.path.exists(mesh_path):
            raise FileNotFoundError(f"Mesh file not found: {mesh_path}")

        vertices = []
        faces = []
        colors = []

        with open(mesh_path, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()

        mode = None
        for line in lines:
            line = line.strip()
            if not line or line.startswith(';'):
                continue

            # Extract Color Names to map the data columns
            if line.startswith('ColorsNames'):
                parts = line.split('=')
                if len(parts) > 1:
                    self.color_names = parts[1].split()
                continue

            if '[VerticesSection]' in line:
                mode = 'v'
                continue
            if '[TrianglesSection]' in line:
                mode = 'f'
                continue
            if '[VerticesColorsSection]' in line:
                mode = 'c'
                continue
            if '[VerticesAttributesSection]' in line:
                mode = 'a'  # Prevent parsing into colors
                continue

            # Parse data lines
            if '=' in line:
                parts = line.split('=', 1)
                if len(parts) == 2:
                    idx_str = parts[0].strip()
                    vals_str = parts[1].strip()
                    try:
                        idx = int(idx_str)
                        vals = vals_str.split()
                        if mode == 'v' and len(vals) >= 3:
                            vertices.append([float(vals[0]), float(vals[1]), float(vals[2])])
                        elif mode == 'f' and len(vals) >= 3:
                            # Carto files index from 1, Python uses 0-based indexing
                            # Reverse winding to fix normal direction (show outside surface)
                            faces.append([int(vals[0]) - 1, int(vals[2]) - 1, int(vals[1]) - 1])
                        elif mode == 'c':
                            colors.append([float(v) for v in vals])
                    except (ValueError, IndexError):
                        continue  # Ignore invalid lines

        self.mesh_vertices = np.array(vertices)
        self.mesh_triangles = np.array(faces)
        self.mesh_colors = np.array(colors) if colors else None
        print(
            f"Mesh loaded: {len(vertices)} vertices, {len(faces)} triangles, {len(self.points)} points, {len(colors)} color entries."
        )

        return self._get_map_data()

    def _points_to_dataframe(self):
        data = []
        for p in self.points:
            pos = p.position
            data.append({
                'ID': p.id,
                'X': pos[0] if pos is not None else None,
                'Y': pos[1] if pos is not None else None,
                'Z': pos[2] if pos is not None else None,
            })
        return pd.DataFrame(data).dropna(subset=['X', 'Y', 'Z'])  # Odrzucamy puste punkty

    def _get_map_data(self):
        """Returns a dictionary with full data (mesh + points)"""
        if self.mesh_vertices is None:
            print("Warning: Mesh not loaded yet. Call load_mesh() first.")

        return {
            'vertices': self.mesh_vertices,
            'triangles': self.mesh_triangles,
            'color_names': self.color_names,
            'colors_mesh': self.mesh_colors,
            'points': self._points_to_dataframe(),
        }

    def modify_model(self, mesh_data, rotate=(0.0, 0.0, 0.0), translate=(0.0, 0.0, 0.0), scale=1.0):
        """Return modified mesh vertices or modified mesh data.

        The method applies scaling, rotation (degrees around X/Y/Z axes), and
        translation to the mesh vertices. If max_distance threshold is set,
        vertices beyond that distance are filtered out.

        Args:
            mesh_data: either a dict returned by load_mesh() containing
                'vertices', or a plain (N, 3) vertices array-like.
            rotate: angles in degrees as (rx, ry, rz).
            translate: offset vector as (dx, dy, dz).
            scale: uniform scale or 3-element scale factors.

        Returns:
            If input is a dict, returns a new dict with transformed 'vertices'
            (and filtered if max_distance is set).
            If input is an array-like vertices list, returns transformed
            vertices as an np.ndarray of shape (N, 3) (filtered if max_distance is set).
        """
        if isinstance(mesh_data, dict):
            if 'vertices' not in mesh_data:
                raise ValueError('mesh_data dict must contain "vertices"')
            vertices = np.asarray(mesh_data['vertices'], dtype=np.float64)
            is_dict_input = True
        else:
            vertices = np.asarray(mesh_data, dtype=np.float64)
            is_dict_input = False

        if vertices.ndim != 2 or vertices.shape[1] != 3:
            raise ValueError('mesh_data must contain an (N, 3) vertices array')

        scale_arr = np.asarray(scale, dtype=np.float64)
        if scale_arr.ndim == 0:
            scale_arr = np.full(3, scale_arr)
        elif scale_arr.shape != (3,):
            raise ValueError('scale must be a scalar or a 3-element vector')

        vertices = vertices * scale_arr

        rotate_arr = np.asarray(rotate, dtype=np.float64)
        if rotate_arr.shape != (3,):
            raise ValueError('rotate must be a 3-element vector of degrees')

        rx, ry, rz = np.deg2rad(rotate_arr)
        cosx, sinx = np.cos(rx), np.sin(rx)
        cosy, siny = np.cos(ry), np.sin(ry)
        cosz, sinz = np.cos(rz), np.sin(rz)

        rx_mat = np.array([[1.0, 0.0, 0.0],
                           [0.0, cosx, -sinx],
                           [0.0, sinx, cosx]], dtype=np.float64)
        ry_mat = np.array([[cosy, 0.0, siny],
                           [0.0, 1.0, 0.0],
                           [-siny, 0.0, cosy]], dtype=np.float64)
        rz_mat = np.array([[cosz, -sinz, 0.0],
                           [sinz, cosz, 0.0],
                           [0.0, 0.0, 1.0]], dtype=np.float64)

        vertices = vertices @ rx_mat.T
        vertices = vertices @ ry_mat.T
        vertices = vertices @ rz_mat.T

        translate_arr = np.asarray(translate, dtype=np.float64)
        if translate_arr.shape != (3,):
            raise ValueError('translate must be a 3-element vector')

        vertices = vertices + translate_arr

        # Apply max_distance filtering if threshold is set
        if self.max_distance is not None and vertices.size > 0:
            dists = np.linalg.norm(vertices, axis=1)
            keep_mask = dists <= float(self.max_distance)
            vertices = vertices[keep_mask]
            
            # If dict input with triangles, filter triangles too
            if is_dict_input and 'triangles' in mesh_data:
                triangles = np.asarray(mesh_data['triangles'], dtype=np.int32)
                keep_indices = np.where(keep_mask)[0]
                idx_map = np.full(mesh_data['vertices'].shape[0], -1, dtype=int)
                idx_map[keep_indices] = np.arange(len(keep_indices), dtype=int)
                
                new_triangles = idx_map[triangles]
                valid_tri_mask = np.all(new_triangles >= 0, axis=1)
                triangles = new_triangles[valid_tri_mask]
                
                modified = mesh_data.copy()
                modified['vertices'] = vertices
                modified['triangles'] = triangles
                return modified

        if is_dict_input:
            modified = mesh_data.copy()
            modified['vertices'] = vertices
            return modified
        return vertices

    def setMaxDistance(self, max_distance):
        """Set the maximum distance threshold for filtering vertices.
        
        When set to a positive value, vertices further than this distance from
        the origin will be filtered out during modify() operations. Pass None
        to disable filtering.
        
        Args:
            max_distance: numeric threshold or None to disable filtering.
        """
        if max_distance is None:
            self.max_distance = None
        else:
            self.max_distance = float(max_distance)

    def get_electrode_types(self):
        """Return available electrode types detected from file_paths."""
        # prefer new hierarchical structure
        electrode_section = self.file_paths.get('Electrode_Positions')
        if isinstance(electrode_section, dict):
            return sorted(electrode_section.keys())

        # Fallback for legacy style
        electrode_types = []
        for key, value in self.file_paths.items():
            if key in ('xml_path', 'mesh_path', 'ECG'):
                continue
            files = value if isinstance(value, (list, tuple)) else [value] if value else []
            if any('Eleclectrode_Positions' in os.path.basename(str(f)) for f in files):
                electrode_types.append(key)
        return electrode_types

    def load_electrodes(self, electrode_type, file_names=None, on_annotation_enable=False):
        """Load electrodes position points for a selected electrode type.

        Args:
            electrode_type: key under self.file_paths.
            file_names: optional list or single entry to override default files.
            on_annotation_enable: if False, ignore files ending with '_OnAnnotation.txt'.
        """

        if electrode_type is None:
            raise ValueError('electrode_type must be specified')

        if file_names is None:
            # prefer new structured data_paths
            electrode_section = self.file_paths.get('Electrode_Positions', {})
            if isinstance(electrode_section, dict) and electrode_type in electrode_section:
                file_names = electrode_section.get(electrode_type, [])
            else:
                files = self.file_paths.get(electrode_type)
                if files is None:
                    legacy_key = f"{electrode_type}_CONNECTOR_Eleclectrode_Positions"
                    files = self.file_paths.get(legacy_key)
                    if files is None:
                        return {}
                file_names = files

        if isinstance(file_names, str):
            file_names = [file_names]

        electrode_data = {}

        for entry in file_names:
            basename = os.path.basename(entry)
            if 'Eleclectrode_Positions' not in basename:
                continue

            if not on_annotation_enable and basename.endswith('_OnAnnotation.txt'):
                continue

            path = entry if os.path.isabs(entry) else os.path.join(self.data_path, entry)
            if not os.path.exists(path):
                continue

            with open(path, 'r', encoding='utf-8', errors='replace') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('Eleclectrode_Positions') or line.startswith('Electrode#'):
                        continue

                    parts = re.split(r'\s+', line)
                    if len(parts) < 5:
                        continue

                    try:
                        eid = int(parts[0])
                        t = float(parts[1])
                        x = float(parts[2])
                        y = float(parts[3])
                        z = float(parts[4])
                    except ValueError:
                        continue

                    key_data = electrode_data.setdefault(eid, [])
                    key_data.append((t, x, y, z))

        return electrode_data

    def add_electrodes(self, file_names=None):
        """Deprecated compatibility wrapper with full scan, not for lazy mode."""
        electrode_type = None
        if isinstance(file_names, str):
            # if passed a named type shortcut
            if file_names in self.file_paths:
                electrode_type = file_names
                file_names = None

        if electrode_type is not None:
            return self.load_electrodes(electrode_type)

        # build generic fallback, matching all electrode types
        electrode_data = {}
        for et in self.get_electrode_types():
            electrode_data.update(self.load_electrodes(et))
        return electrode_data

