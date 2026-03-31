import os

class CartoMapDataPaths:
    def __init__(self, data_path, measurement_name):
        self.data_path = os.path.abspath(data_path)
        self.measurement_name = measurement_name

        self.xml_path = os.path.join(self.data_path, f"{measurement_name}_Points_Export.xml")
        self.mesh_path = os.path.join(self.data_path, f"{measurement_name}.mesh")

        self.file_paths = {
            'xml_path': self.xml_path,
            'mesh_path': self.mesh_path,
        }

        self._scan_position_files()

        self.master_xml = self.xml_path
        self.base_dir = self.data_path

    def _scan_position_files(self):
        """Scan data_path for "{measurement_name}_{some_name}_{some_numbers}.txt" files."""
        import re

        if not os.path.isdir(self.data_path):
            return

        # Keep names like measurement_name_a_b_c_123.txt as key='a_b_c'
        matcher = re.compile(rf"^{re.escape(self.measurement_name)}_(.+?)_(\d+)\.txt$")
        grouped = {}

        for fn in os.listdir(self.data_path):
            if not fn.lower().endswith('.txt'):
                continue
            m = matcher.match(fn)
            if not m:
                continue

            key = m.group(1)
            grouped.setdefault(key, []).append(fn)

        # sort lists for stable order
        for key, values in grouped.items():
            grouped[key] = sorted(values)

        self.file_paths.update(grouped)

