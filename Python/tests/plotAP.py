import os, sys
import numpy as np
import plotly.graph_objects as go

from pathlib import Path
from dotenv import load_dotenv
load_dotenv()
data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()

def plot_ap(ap_signals, signal_names=None, title="AP signals"):
    if signal_names is None:
        signal_names = [f"Signal {i+1}" for i in range(len(ap_signals))]
    
    fig = go.Figure()
    for i, ap in enumerate(ap_signals):
        fig.add_trace(go.Scatter(x=list(range(len(ap))), y=ap, mode='lines', name=signal_names[i]))
    
    fig.update_layout(
        title=title,
        xaxis_title='ms',
        yaxis_title='mV'
    )
    fig.show()

signal_names=["Endocardial", "Epicardial", "Midmyocardial"]
signals = [np.loadtxt(os.path.join(data_path, 'AP', f'AP_{name.lower()}.txt'))[32000:] for name in signal_names]

plot_ap(signals, signal_names, title="AP signals")