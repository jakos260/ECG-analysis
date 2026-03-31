import os, sys
import json
from pathlib import Path
from dotenv import load_dotenv

def dict2json(d, filename):
    with open(filename, 'w') as f:
        json.dump(d, f, indent=4)

load_dotenv()
root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))
from carto3reader.CartoMap import CartoMap

data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
sample_data_path = os.path.join(data_path, 'Carto3Data', 'Patient 2025_10_20', 'VT1', 'Export_VT1-02_02_2026-12-24-27')

measurement_name = '1-VT'

carto_map = CartoMap(sample_data_path, measurement_name)
electrodes_types = carto_map.get_electrode_types()
print(electrodes_types)

cs = carto_map.load_electrodes(electrodes_types[0])
print(cs.keys())
# print(cs.keys())
# dict2json(carto_map.file_paths, f'readCartoData/_data/{measurement_name}_file_paths.json')
# dict2json(electrodes, f'readCartoData/_data/{measurement_name}_electrodes.json')

