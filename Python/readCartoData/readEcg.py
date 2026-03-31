import os, sys
import numpy as np
import xml.etree.ElementTree as ET
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

from pathlib import Path
from dotenv import load_dotenv

root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))                         # add helpers
from carto3reader.CartoMap import CartoMap

load_dotenv()

def plotEcg(carto_map):
    try:
        ecg_fig = carto_map.ecg_reader.plot()
        return ecg_fig
    except Exception as e:
        print("Error reading ECG data:", e)
        return None

if __name__ == "__main__":
    
    data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()
    sample_data_path = os.path.join(data_path, 'Carto3Data', 'Patient 2025_10_20', 'VT1', 'Export_VT1-02_02_2026-12-24-27')

    measurement_name = '1-VT' # 114 points
    # measurement_name = '2-AT LV'
    # measurement_name = '2-1-ReAT LV'
    # measurement_name = '1-1-1-PATTERN-STYMU'
    # measurement_name = '1-1-PATTERN BEZ STYMU'
    carto_map = CartoMap(sample_data_path, measurement_name)
    plotEcg(carto_map)







