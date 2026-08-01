import os, sys
import traceback
import numpy as np
from scipy.signal import find_peaks
import json
from pathlib import Path

root_path = Path(__file__).resolve().parent.parent
sys.path.append(str(root_path))
from helpers.helpers import loadmat

# organising function
def process_ecg_from_patient_folder(patient_name, dataset_path, subfolder='signals'):
    """
    Process ECG signals from a specified subfolder within the output directory.
    
    Parameters:
    - patient_name: The name of the patient folder.
    - dataset_path: The path to the dataset directory.
    - subfolder: The specific subfolder containing ECG signal files (e.g., 'signals').
    
    This function iterates through each patient folder, processes the ECG signals,
    and saves the processed results back to the respective patient folder.
    """        

    signal_path = os.path.join(dataset_path, subfolder)
    if not os.path.isdir(signal_path):
        print(f"  Skipping {patient_name}: no '{subfolder}' directory found at {signal_path}")
        return

    patient_info = {}
    for signal_file in sorted(os.listdir(signal_path)):
        if not signal_file.endswith('.bsm'):
            continue
        full_signal_path = os.path.join(signal_path, signal_file)
        try:
            bsm, _ = loadmat(full_signal_path)
        except Exception as e:
            print(f"  Error loading {full_signal_path}: {e}")
            print(traceback.format_exc())
            continue

        print(f"Processing {full_signal_path} with shape: {bsm.shape}")
        sampling_rate = 1000
        try:
            processed_signals, valid_beats, beats_begin = process_ecg(bsm, sampling_rate)
        except Exception as e:
            print(f"  Error processing {signal_file}: {e}")
            print(traceback.format_exc())
            continue

        signal_name = os.path.splitext(signal_file)[0]

        if isinstance(valid_beats, np.ndarray):
            beats_list = valid_beats.tolist()
        else:
            beats_list = list(valid_beats)
        patient_info[signal_name] = {
                        "beats": beats_list
                    }

        json_file_path = os.path.join(signal_path, "info.json")
        with open(json_file_path, 'w', encoding='utf-8') as f:
            json.dump(patient_info, f, indent=4, ensure_ascii=False)

        print(f"  Saved info.json for {patient_name} in {signal_path}")
    else:
        print(f"  No .bsm files found in {signal_path}; skipped info.json")

def process_ecg(signal, sampling_rate):

    processed_signals = np.zeros_like(signal)
    all_peaks = []
    for i in range(signal.shape[0]):
        tmp = signal[i,:]
        tmp = tmp - np.mean(tmp)  # Remove DC offset
        tmp = tmp / np.max(np.abs(tmp))  # Normalize
        [are_peaks_positive, positive_peaks, negative_peaks] = check_if_peaks_positive(tmp, sampling_rate)
        peaks = positive_peaks if are_peaks_positive else negative_peaks

        processed_signals[i,:] = tmp
        all_peaks.append(peaks)
    
    valid_beats = find_consensus_beats(
        np.concatenate(all_peaks),
        tolerance_samples=int(0.05*sampling_rate),
        min_votes=10)
    
    beats_begin = valid_beats - int(0.3*sampling_rate)
    return processed_signals, valid_beats, beats_begin

def check_if_peaks_positive(signal, sampling_rate):
    positive_peaks = find_peaks_wrapper(signal, sampling_rate)
    negative_peaks = find_peaks_wrapper(-signal, sampling_rate)
    if len(positive_peaks) == 0 and len(negative_peaks) == 0:
        return None, [], []
    if len(positive_peaks) == 0:
        return False, [], negative_peaks
    if len(negative_peaks) == 0:
        return True, positive_peaks, []
    are_peaks_positive = np.std(np.diff(positive_peaks)) < np.std(np.diff(negative_peaks))
    # print(f"Positive peaks std: {np.diff(positive_peaks)}, Negative peaks std: {np.diff(negative_peaks)}")
    return are_peaks_positive, positive_peaks, negative_peaks

def find_peaks_wrapper(signal, sampling_rate):
    peaks, _ = find_peaks(
        signal,
        height=0.4,
        distance=sampling_rate*0.6,
        prominence=0.25
        )
    return peaks

def find_consensus_beats(all_peaks, tolerance_samples, min_votes, min_distance_samples=500):
    """
    Znajduje docelowe pozycje beatów z globalną weryfikacją rytmu.
    
    Parametry:
    - min_distance_samples: minimalny fizjologiczny odstęp między pobudzeniami (np. 250 próbek dla 1000 Hz = 250 ms)
    """
    # --- ETAP 1: Klastrowanie (grupowanie bliskich detekcji) ---
    sorted_peaks = np.sort(all_peaks)
    
    if len(sorted_peaks) == 0:
        return np.array([], dtype=int)

    clusters = []
    current_cluster = [sorted_peaks[0]]
    
    for peak in sorted_peaks[1:]:
        if peak - current_cluster[-1] <= tolerance_samples:
            current_cluster.append(peak)
        else:
            clusters.append(current_cluster)
            current_cluster = [peak]
    clusters.append(current_cluster)
    
    # --- ETAP 2: Wybór kandydatów ---
    candidates = []
    votes = []
    
    for cluster in clusters:
        if len(cluster) >= min_votes:
            # POPRAWKA: Rzutowanie mediany na int, aby indeksy były poprawne
            candidates.append(int(np.median(cluster)))
            votes.append(len(cluster)) 
            
    candidates = np.array(candidates)
    votes = np.array(votes)
    
    if len(candidates) < 2:
        return candidates

    # --- ETAP 3: Globalne rozwiązywanie konfliktów (Iteracyjna eliminacja) ---
    # Zamiast iść od lewej do prawej, szukamy globalnie najciaśniejszych par
    # i usuwamy z nich słabszego kandydata, aż wszystkie odstępy będą >= min_distance_samples
    
    while True:
        rr_intervals = np.diff(candidates)
        
        # Jeśli nie ma już konfliktów (wszystkie dystanse są OK), przerywamy pętlę
        if len(rr_intervals) == 0 or np.min(rr_intervals) >= min_distance_samples:
            break
            
        # Znajdujemy indeks pierwszego konfliktu (najmniejszy dystans)
        conflict_idx = np.argmin(rr_intervals)
        
        # Porównujemy "głosy" dwóch klastrów biorących udział w konflikcie
        # conflict_idx to indeks pierwszego z pary, conflict_idx + 1 to indeks drugiego
        if votes[conflict_idx] < votes[conflict_idx + 1]:
            loser_idx = conflict_idx
        else:
            loser_idx = conflict_idx + 1
            
        # Usuwamy przegranego kandydata z obu tablic
        candidates = np.delete(candidates, loser_idx)
        votes = np.delete(votes, loser_idx)
            
    return candidates


if __name__ == "__main__":
    from dotenv import load_dotenv
    import plotly.graph_objects as go
    load_dotenv()

    root_path = Path(__file__).resolve().parent.parent

    data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()

    patient_name, measurement_name = 'IKEM_Pat001', '3_2019-11-15-12-36-17'
    # patient_name, measurement_name = 'IKEM_Pat005', '7_AVd_150_VVd_-40_2019-12-9-9-50-21'
    signal_path = os.path.join(data_path, 'Dataset', patient_name, 'signals', f'{measurement_name}.bsm')

    
    bsm, _ = loadmat(signal_path)
    print(f"Loaded ECG signal shape: {bsm.shape}")

    sampling_rate = 1000
    processed_signals, valid_beats, beats_begin = process_ecg(bsm, sampling_rate)
    print(f"Detected {valid_beats.shape[0]} valid beats at positions: {valid_beats}")
    
    fig = go.Figure()
    for i in range(bsm.shape[0]):
        # fig.add_trace(go.Scatter(x=list(range(len(bsm[i,:]))), y=bsm[i,:], name=f'Lead {i}'))
        fig.add_trace(go.Scatter(x=np.arange(processed_signals.shape[1]//1000*sampling_rate), y=processed_signals[i,:]-i, name=f'Processed {i}'))
        fig.add_trace(go.Scatter(x=valid_beats, y=processed_signals[i,valid_beats]-i, marker=dict(color='black', size=8), mode='markers'))
        fig.add_trace(go.Scatter(x=beats_begin, y=processed_signals[i,beats_begin]-i, marker=dict(color='blue', size=8), mode='markers'))

    fig.update_xaxes(title_text="[ms]")
    fig.update_yaxes(title_text="Amplitude")
    fig.update_layout(title=f'signal {signal_path}')
    fig.show()