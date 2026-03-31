import os
import types

class CartoMapDataPaths:
    def __init__(self, data_path, measurement_name):
        self.data_path = os.path.abspath(data_path)
        self.measurement_name = measurement_name

        self.xml_path = os.path.join(self.data_path, f"{measurement_name}_Points_Export.xml")
        self.mesh_path = os.path.join(self.data_path, f"{measurement_name}.mesh")

        self.file_paths = {
            'xml_path': self.xml_path,
            'mesh_path': self.mesh_path,
            'Electrode_Positions': {},
            'Sensor_Positions': {},
            'ECG': [],
        }

        self._scan_position_files()
        self._update_data_paths_namespace()

        self.master_xml = self.xml_path
        self.base_dir = self.data_path

    def _update_data_paths_namespace(self):
        def _to_namespace(value):
            if isinstance(value, dict):
                ns = types.SimpleNamespace()
                for k, v in value.items():
                    setattr(ns, k, _to_namespace(v))
                return ns
            return value

        positions = {
            'Electrode': self.file_paths.get('Electrode_Positions', {}),
            'Sensor': self.file_paths.get('Sensor_Positions', {}),
        }

        self.data_paths = types.SimpleNamespace(
            Positions=_to_namespace(positions),
            Electrode_Positions=_to_namespace(self.file_paths.get('Electrode_Positions', {})),
            Sensor_Positions=_to_namespace(self.file_paths.get('Sensor_Positions', {})),
            ECG=self.file_paths.get('ECG', []),
            xml_path=self.xml_path,
            mesh_path=self.mesh_path,
        )

    def _scan_position_files(self):
        """Scan data_path for "{measurement_name}_{some_name}_{some_numbers}.txt" files."""
        import re

        if not os.path.isdir(self.data_path):
            return

        matcher = re.compile(rf"^{re.escape(self.measurement_name)}_(.+?)_(\d+)\.txt$")

        for fn in os.listdir(self.data_path):
            if not fn.lower().endswith('.txt'):
                continue

            m = matcher.match(fn)
            if not m:
                continue

            key = m.group(1)

            if key.upper().startswith('ECG'):
                if fn not in self.file_paths['ECG']:
                    self.file_paths['ECG'].append(fn)
                continue

            # Normalize category names and device names
            # Example key patterns:
            # CS_CONNECTOR_Eleclectrode_Positions
            # MCC_DX_CONNECTOR_Sensor_Positions
            # QUAD_A_CONNECTOR_Eleclectrode_Positions_OnAnnotation
            m2 = re.match(r'^(?P<device>.+?)_CONNECTOR_(?P<kind>Eleclectrode|Sensor)_Positions(?:_OnAnnotation)?$', key)
            if not m2:
                continue

            device = m2.group('device').strip('_')
            kind = m2.group('kind')

            if kind == 'Eleclectrode':
                category = 'Electrode_Positions'
            else:
                category = 'Sensor_Positions'

            target = self.file_paths.setdefault(category, {})
            target.setdefault(device, []).append(fn)

        # sort lists for stable order
        for dev_group in ('Electrode_Positions', 'Sensor_Positions'):
            for device_name, files in self.file_paths.get(dev_group, {}).items():
                files.sort()
        self.file_paths['ECG'].sort()


