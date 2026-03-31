import os
import re
import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go


class EcgReader:
    def __init__(self, carto_map):
        self.carto_map = carto_map
        self.figure = None

    def get_files(self):
        """Return list of ECG export files from CartoMap.file_paths."""
        files = self.carto_map.file_paths.get('ECG', [])
        if not files:
            files = self.carto_map.file_paths.get('ECG_Export', [])

        if isinstance(files, str):
            files = [files]
        return list(files) if files is not None else []

    def _resolve_path(self, filename):
        if os.path.isabs(filename):
            return filename
        return os.path.join(self.carto_map.data_path, filename)

    def read_data(self, file_index=0):
        """Read ECG data and return pandas.DataFrame."""
        ecg_files = self.get_files()
        if not ecg_files:
            raise FileNotFoundError("No ECG_Export files found in carto_map.file_paths")

        if file_index < 0 or file_index >= len(ecg_files):
            raise IndexError(f"ECG file index out of range: {file_index}")

        filename = ecg_files[file_index]
        filepath = self._resolve_path(filename)

        if not os.path.exists(filepath):
            raise FileNotFoundError(f"ECG file not found: {filepath}")

        headers = None
        rows = []

        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('ECG_Export') or line.startswith('Raw ECG'):
                    continue

                parts = re.split(r'\s+', line)
                is_numeric = all(re.fullmatch(r'-?\d+(\.\d+)?', p) for p in parts if p != '')

                if not is_numeric:
                    headers = parts
                    continue

                row = [float(p) for p in parts if p != '']
                rows.append(row)

        if not rows:
            raise ValueError(f"No numeric ECG data found in file: {filepath}")

        df = pd.DataFrame(rows)

        if headers and len(headers) == df.shape[1]:
            df.columns = headers

        return df

    def plot(self, file_index=0, channel=None, max_channels=8):
        """Plot ECG (plotly) from ECG export data."""
        df = self.read_data(file_index=file_index)

        if channel is None:
            columns = df.columns if len(df.columns) <= max_channels else df.columns[:max_channels]
        elif isinstance(channel, int):
            columns = [df.columns[channel]]
        else:
            columns = [channel]

        fig = go.Figure()
        for col in columns:
            fig.add_trace(go.Scatter(x=df.index, y=df[col], mode='lines', name=str(col)))

        fig.update_layout(
            title=f"ECG data: {self.get_files()[file_index]}",
            xaxis_title='Sample',
            yaxis_title='Value',
            height=600,
        )

        self.figure = fig
        fig.show()
        return fig

    def close(self):
        """Close the current plotted figure if available."""
        print("Closing browser tab doesn't work")
        if self.figure is None:
            return

        try:
            if hasattr(pio, 'close'):
                pio.close(self.figure)
            elif hasattr(self.figure, 'close'):
                self.figure.close()
        except Exception:
            pass

        self.figure = None
