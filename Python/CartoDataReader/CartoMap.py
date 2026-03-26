import os
import pandas as pd
import xml.etree.ElementTree as ET
import numpy as np

class CartoPoint:
    def __init__(self, point_id, xml_path, base_dir):
        self.id = point_id
        self.xml_path = xml_path
        self.base_dir = base_dir
        
        # Podstawowe parametry (wczytywane od razu z XML punktu)
        self.bipolar = None
        self.unipolar = None
        self.ecg_filename = None
        self.pos_filename = None
        
        # Dane "ciężkie" (ładowane na żądanie)
        self._position_df = None
        self._ecg_df = None
        
        self._parse_point_xml()

    def _parse_point_xml(self):
        """Parsuje plik XML konkretnego punktu (Pxxx.xml)"""
        tree = ET.parse(self.xml_path)
        root = tree.getroot()
        
        # Pobieranie woltażu
        v_node = root.find('Voltages')
        if v_node is not None:
            self.bipolar = float(v_node.get('Bipolar', 0))
            self.unipolar = float(v_node.get('Unipolar', 0))
            
        # Szukanie plików powiązanych
        ecg_node = root.find('ECG')
        if ecg_node is not None:
            self.ecg_filename = ecg_node.get('FileName')
            
        for conn in root.findall('.//Connector'):
            for attr_name, attr_val in conn.attrib.items():
                if 'OnAnnotation' in attr_val:
                    self.pos_filename = attr_val
                    break

    def _read_carto_tsv(self, path, **kwargs):
        """Read a Carto-exported TSV file.

        Some Carto files include a single title line (e.g. "Sensor_Positions_2.0")
        before the tab-separated header row. In that case, pandas must use
        header=1 (skip the title line) rather than skiprows=2 (which would
        treat the first data row as the header and result in an empty dataframe).
        """
        # Detect whether the first line is a title line (no tabs) and the second line is the true header.
        header = 0
        try:
            with open(path, 'r', encoding='utf-8', errors='replace') as f:
                first = f.readline()
                second = f.readline()
        except FileNotFoundError:
            raise

        if first and '\t' not in first and second and '\t' in second:
            header = 1

        # Ensure the first column is treated as data (not an index). Some Carto
        # exports include a leading numeric column that pandas can mistakenly
        # treat as the DataFrame index if index_col is left to pandas' default.
        return pd.read_csv(path, sep='\t', header=header, index_col=False, **kwargs)

    @property
    def position(self):
        """Zwraca wektor [X, Y, Z] dla momentu annotacji (Lazy Loading)"""
        if self._position_df is None and self.pos_filename:
            path = os.path.join(self.base_dir, self.pos_filename)
            df = self._read_carto_tsv(path, usecols=['X', 'Y', 'Z'], nrows=1)
            if not df.empty:
                self._position_df = df.iloc[0].values
        return self._position_df

    def ecg(self):
        """Wczytuje i zwraca pełny sygnał EKG dla punktu"""
        if self._ecg_df is None and self.ecg_filename:
            path = os.path.join(self.base_dir, self.ecg_filename)
            self._ecg_df = self._read_carto_tsv(path)
        return self._ecg_df

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
        
        self._parse_point_xml()

    def _parse_point_xml(self):
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

    def _read_carto_tsv(self, path, **kwargs):
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
                df = self._read_carto_tsv(path, usecols=['X', 'Y', 'Z'], nrows=1)
                if not df.empty:
                    self._position_df = df.iloc[0].values
        return self._position_df

    def ecg(self):
        if self._ecg_df is None and self.ecg_filename:
            path = os.path.join(self.base_dir, self.ecg_filename)
            if os.path.exists(path):
                self._ecg_df = self._read_carto_tsv(path)
        return self._ecg_df

class CartoMap:
    def __init__(self, master_xml_path):
        self.master_xml = master_xml_path
        self.base_dir = os.path.dirname(master_xml_path)
        self.points = []
        
        # Nowe atrybuty na siatkę 3D
        self.mesh_vertices = None
        self.mesh_triangles = None
        
        self._load_master()

    def _load_master(self):
        tree = ET.parse(self.master_xml)
        for p_entry in tree.getroot().findall('Point'):
            p_id = p_entry.get('ID')
            p_xml = os.path.join(self.base_dir, p_entry.get('File_Name'))
            
            if os.path.exists(p_xml):
                self.points.append(CartoPoint(p_id, p_xml, self.base_dir))

    def load_mesh(self, mesh_filename):
        """Loads heart anatomy from .mesh file"""
        mesh_path = os.path.join(self.base_dir, mesh_filename)
        if not os.path.exists(mesh_path):
            raise FileNotFoundError(f"Mesh file not found: {mesh_path}")
        
        vertices = []
        faces = []
        
        with open(mesh_path, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
            
        mode = None
        for line in lines:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
                
            if '[VerticesSection]' in line:
                mode = 'v'
                continue
            if '[TrianglesSection]' in line:
                mode = 'f'
                continue
            
            parts = line.split()
            if len(parts) >= 5 and '=' in parts:  # Ensure it's a data line with index = values
                try:
                    if mode == 'v':
                        vertices.append([float(parts[2]), float(parts[3]), float(parts[4])])
                    elif mode == 'f':
                        # Carto files index from 1, Python uses 0-based indexing
                        # Reverse winding to fix normal direction (show outside surface)
                        faces.append([int(parts[2])-1, int(parts[4])-1, int(parts[3])-1])
                except (ValueError, IndexError):
                    continue  # Ignore invalid lines

        # Convert to NumPy arrays
        self.mesh_vertices = np.array(vertices)
        self.mesh_triangles = np.array(faces)
        print(f"Mesh loaded: {len(vertices)} vertices, {len(faces)} triangles.")
        return self.mesh_vertices, self.mesh_triangles

    def to_dataframe(self):
        data = []
        for p in self.points:
            pos = p.position
            data.append({
                'ID': p.id,
                'X': pos[0] if pos is not None else None,
                'Y': pos[1] if pos is not None else None,
                'Z': pos[2] if pos is not None else None,
                'Bipolar': p.bipolar,
                'Unipolar': p.unipolar
            })
        return pd.DataFrame(data).dropna(subset=['X', 'Y', 'Z']) # Odrzucamy puste punkty

    def get_map_data(self):
        """Returns a dictionary with full data (mesh + points)"""
        if self.mesh_vertices is None:
            print("Warning: Mesh not loaded yet. Call load_mesh() first.")
            
        return {
            'vertices': self.mesh_vertices,
            'triangles': self.mesh_triangles,
            'points_df': self.to_dataframe()
        }