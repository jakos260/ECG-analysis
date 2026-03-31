import os, sys
import json
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()
root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))
from carto3reader.CartoMap import CartoMap

data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
sample_data_path = os.path.join(data_path, 'Carto3Data', 'Patient 2025_10_20', 'VT1', 'Export_VT1-02_02_2026-12-24-27')

# measurement_name = '1-VT'
measurement_name = '2-AT LV'

carto_map = CartoMap(sample_data_path, measurement_name)
with open(f'readCartoData/_data/{measurement_name}_file_paths.json', 'w') as f:
    json.dump(carto_map.file_paths, f, indent=4)