import os
import pandas as pd
import xml.etree.ElementTree as ET
import numpy as np

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

class CartoMap:
    def __init__(self, master_xml_path):
        self.master_xml = master_xml_path
        self.base_dir = os.path.dirname(master_xml_path)
        self.points = []
        
        self.mesh_vertices = None
        self.mesh_triangles = None
        
        # New attributes for vertex colors
        self.mesh_colors = None
        self.color_names = []
        
        self.__load_master()

    def __load_master(self):
        tree = ET.parse(self.master_xml)
        for p_entry in tree.getroot().findall('Point'):
            p_id = p_entry.get('ID')
            p_xml = os.path.join(self.base_dir, p_entry.get('File_Name'))
            
            if os.path.exists(p_xml):
                self.points.append(CartoPoint(p_id, p_xml, self.base_dir))

    def load_mesh(self, mesh_filename):
        """Loads heart anatomy and associated color data from .mesh file"""
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
                mode = 'a' # Prevent parsing into colors
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
                            faces.append([int(vals[0])-1, int(vals[2])-1, int(vals[1])-1])
                        elif mode == 'c':
                            colors.append([float(v) for v in vals])
                    except (ValueError, IndexError):
                        continue  # Ignore invalid lines

        self.mesh_vertices = np.array(vertices)
        self.mesh_triangles = np.array(faces)
        self.mesh_colors = np.array(colors) if colors else None
        print(f"Mesh loaded: {len(vertices)} vertices, {len(faces)} triangles, {len(self.points)} points, {len(colors)} color entries.")
        
        return self.__get_map_data()

    def __points_to_dataframe(self):
        data = []
        for p in self.points:
            pos = p.position
            data.append({
                'ID': p.id,
                'X': pos[0] if pos is not None else None,
                'Y': pos[1] if pos is not None else None,
                'Z': pos[2] if pos is not None else None
            })
        return pd.DataFrame(data).dropna(subset=['X', 'Y', 'Z']) # Odrzucamy puste punkty

    def __get_map_data(self):
        """Returns a dictionary with full data (mesh + points)"""
        if self.mesh_vertices is None:
            print("Warning: Mesh not loaded yet. Call load_mesh() first.")
            
        return {
            'vertices': self.mesh_vertices,
            'triangles': self.mesh_triangles,
            'color_names': self.color_names,
            'colors_mesh': self.mesh_colors,
            'points': self.__points_to_dataframe()
        }