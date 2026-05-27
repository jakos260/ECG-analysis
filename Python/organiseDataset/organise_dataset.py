import os
import shutil
import json
import re

# Definiowanie ścieżek
ROOT_DIR = r"C:\Users\Admin\Documents\Projects\ecg project\Scripts\data"
BASE_DIR = os.path.join(ROOT_DIR, "bin")
MODELS_DIR = os.path.join(BASE_DIR, "Models")
MAPPER_DIR = os.path.join(BASE_DIR, "Mapper")
ECG_DATA_DIR = os.path.join(MAPPER_DIR, "ECG_DATA")

# Folder docelowy na uporządkowane dane
OUTPUT_DIR = os.path.join(ROOT_DIR, "Dataset")

def get_patient_data(dest_patient_dir, clean_patient_name):
    """
    Buduje strukturę pacjenta, z różną logiką pomiarów dla ECGSIM i IKEM.
    """
    dest_model_dir = os.path.join(dest_patient_dir, "model")
    dest_ecgs_dir = os.path.join(dest_patient_dir, "ecgs")
    dest_beats_dir = os.path.join(dest_patient_dir, "ventricular_beats")
    
    is_ecgsim = "ECGSIM" in clean_patient_name.upper()
    is_ikem = "IKEM" in clean_patient_name.upper()

    # --- CZĘŚĆ: MODEL ---
    ct_path = ""
    ama_bsm_path = ""
    ama_s12_path = ""
    
    if os.path.exists(dest_model_dir):
        for f in os.listdir(dest_model_dir):
            if f == "ventricle.typ":
                ct_path = f
            elif "ventricles2BSM" in f and f.endswith(".mat"):
                ama_bsm_path = f
            elif f == "ventricles2standard_12.mat":
                ama_s12_path = f

    model_dict = {
        "VER_ITRI": "ventricle.tri",
        "DIST": "ventricle.dst3d",
        "DISTsurf": "ventricle.dst2d",
        "DIST2W": "ventricle.dstanis",
        "ADJ2W": "ventricle.adjanis",
        "neigh": "ventricle.adjneigh",
        "ADJ": "ventricle.adj3d",
        "CT": ct_path
    }

    # --- CZĘŚĆ: SIGNALS ---
    
    # Określanie liczby elektrod z pliku refECG zawierającego 'BSM'
    bsm_leads = -1
    if os.path.exists(dest_ecgs_dir):
        for f in os.listdir(dest_ecgs_dir):
            if "BSM" in f and f.endswith(".refECG"):
                full_path = os.path.join(dest_ecgs_dir, f)
                try:
                    with open(full_path, 'r', encoding='utf-8') as file:
                        bsm_leads = sum(1 for line in file if line.strip())
                except Exception:
                    pass
                break

    bsm_leads_loc = {
        "leads_detected": bsm_leads,
        "AMA": ama_bsm_path,
        "elec": "",
        "iecg": "",
        "imap": "",
        "imapLog": ""
    }

    s12_leads_loc = {
        "AMA": ama_s12_path,
        "elec": "",
        "iecg": "",
        "imap": "",
        "imapLog": ""
    }

    # Budowa sekcji POMIARÓW z podziałem na rodzaj pacjenta
    measurements = []
    
    if is_ecgsim:
        bsm_ref = ""
        ecg_ref = ""
        
        # Szukamy plików refECG
        if os.path.exists(dest_ecgs_dir):
            for f in os.listdir(dest_ecgs_dir):
                if "BSM" in f and f.endswith(".refECG"):
                    bsm_ref = f
                elif f == "standard_12.refECG":
                    ecg_ref = f
        
        measurements.append({
            "name": "ecgs",
            "bsm": bsm_ref,
            "bsm_medianecg": "",
            "ecg": ecg_ref
        })

        # Szukamy powiązanych plików .elec (elektrod)
        if bsm_ref:
            bsm_elec = bsm_ref.replace(".refECG", ".elec")
            if os.path.exists(os.path.join(dest_ecgs_dir, bsm_elec)):
                bsm_leads_loc["elec"] = bsm_elec
        if ecg_ref:
            ecg_elec = ecg_ref.replace(".refECG", ".elec")
            if os.path.exists(os.path.join(dest_ecgs_dir, ecg_elec)):
                s12_leads_loc["elec"] = ecg_elec

    elif is_ikem:
        grouped_measurements = {}
        if os.path.exists(dest_ecgs_dir):
            for f in os.listdir(dest_ecgs_dir):
                # Ekstrakcja unikalnej nazwy (bez rozszerzeń i prefiksów)
                if f.endswith(".bsm.medianecg"):
                    m_name = f[:-14]
                    if m_name not in grouped_measurements:
                        grouped_measurements[m_name] = {"bsm": "", "bsm_medianecg": "", "ecg": ""}
                    grouped_measurements[m_name]["bsm_medianecg"] = f
                
                elif f.endswith(".bsm"):
                    m_name = f[:-4]
                    if m_name not in grouped_measurements:
                        grouped_measurements[m_name] = {"bsm": "", "bsm_medianecg": "", "ecg": ""}
                    grouped_measurements[m_name]["bsm"] = f
                    
                elif f.endswith(".ecg"):
                    m_name = f[:-4]
                    if m_name not in grouped_measurements:
                        grouped_measurements[m_name] = {"bsm": "", "bsm_medianecg": "", "ecg": ""}
                    grouped_measurements[m_name]["ecg"] = f
        
        # Przekonwertowanie ze słownika na listę wymaganą w JSON
        for m_name in sorted(grouped_measurements.keys()):
            measurements.append({
                "name": m_name,
                "bsm": grouped_measurements[m_name]["bsm"],
                "bsm_medianecg": grouped_measurements[m_name]["bsm_medianecg"],
                "ecg": grouped_measurements[m_name]["ecg"]
            })

    # Ground truth (sprawdzanie konkretnie user.dep i user.rep)
    dep_path, rep_path = "", ""
    if os.path.exists(dest_beats_dir):
        beat1_dir = os.path.join(dest_beats_dir, "beat1")
        if os.path.exists(beat1_dir):
            for f in os.listdir(beat1_dir):
                if f == "user.dep":
                    dep_path = f
                elif f == "user.rep":
                    rep_path = f

    signals_dict = {
        "bsm_leads_loc": bsm_leads_loc,
        "s12_leads_loc": s12_leads_loc,
        "measurements": measurements,
        "true_dep": dep_path,
        "true_rep": rep_path
    }

    return {
        clean_patient_name: {
            "model": model_dict,
            "signals": signals_dict
        }
    }


def organize_dataset():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    global_geom_info = {
        "patients": []
    }

    # Przejście przez wszystkie foldery w Dataset/Models
    for patient_folder in os.listdir(MODELS_DIR):
        model_path = os.path.join(MODELS_DIR, patient_folder)
        if not os.path.isdir(model_path):
            continue

        clean_patient_name = patient_folder
        if clean_patient_name.endswith("_model"):
            clean_patient_name = clean_patient_name[:-6]

        print(f"Przetwarzam: {patient_folder} -> {clean_patient_name}")

        core_id_match = re.search(r'(Pat\d{3})', patient_folder)
        core_id = core_id_match.group(1) if core_id_match else None

        dest_patient_dir = os.path.join(OUTPUT_DIR, clean_patient_name)
        dest_model_dir = os.path.join(dest_patient_dir, "model")
        dest_ecgs_dir = os.path.join(dest_patient_dir, "ecgs")
        dest_ventricular_dir = os.path.join(dest_patient_dir, "ventricular_beats")

        os.makedirs(dest_model_dir, exist_ok=True)
        os.makedirs(dest_ecgs_dir, exist_ok=True)

        # 1. Kopiowanie z zachowaniem struktury
        for item in os.listdir(model_path):
            src_item = os.path.join(model_path, item)
            if os.path.isdir(src_item):
                if item == "ventricular_beats":
                    os.makedirs(dest_ventricular_dir, exist_ok=True)
                    shutil.copytree(src_item, dest_ventricular_dir, dirs_exist_ok=True)
                elif item == "ecgs":
                    shutil.copytree(src_item, dest_ecgs_dir, dirs_exist_ok=True)
                else:
                    shutil.copytree(src_item, os.path.join(dest_model_dir, item), dirs_exist_ok=True)
            else:
                shutil.copy2(src_item, dest_model_dir)

        # 2. Szukanie plików EKG w Mapper
        if core_id:
            if os.path.exists(MAPPER_DIR):
                for file in os.listdir(MAPPER_DIR):
                    if core_id in file and os.path.isfile(os.path.join(MAPPER_DIR, file)):
                        shutil.copy2(os.path.join(MAPPER_DIR, file), dest_ecgs_dir)

            if os.path.exists(ECG_DATA_DIR):
                for file in os.listdir(ECG_DATA_DIR):
                    if core_id in file and os.path.isfile(os.path.join(ECG_DATA_DIR, file)):
                        file_path = os.path.join(ECG_DATA_DIR, file)
                        shutil.copy2(file_path, dest_ecgs_dir)

        # Tworzenie info.json
        info_data = {"qrs_begin": -1, "t_begin": -1, "BSM_leads_number": -1}
        with open(os.path.join(dest_patient_dir, "info.json"), "w", encoding='utf-8') as f:
            json.dump(info_data, f, indent=4)

        # 3. CZYSZCZENIE NAZW PLIKÓW (z prefiksów)
        prefixes_to_remove = [
            f"{patient_folder}_",      
            f"{clean_patient_name}_",  
            f"{core_id}_"              
        ]
        prefixes_to_remove = [p for p in prefixes_to_remove if p != "_"]

        for root, dirs, files in os.walk(dest_patient_dir):
            for file in files:
                new_name = file
                for prefix in prefixes_to_remove:
                    if new_name.startswith(prefix):
                        new_name = new_name[len(prefix):].lstrip(" _")
                        break 
                
                if new_name != file:
                    old_path = os.path.join(root, file)
                    new_path = os.path.join(root, new_name)
                    if not os.path.exists(new_path):
                        os.rename(old_path, new_path)

        # 4. Generowanie słownika i dodawanie do głównego JSON
        patient_dictionary = get_patient_data(dest_patient_dir, clean_patient_name)
        global_geom_info["patients"].append(patient_dictionary)

    # 5. Zapis
    global_info_path = os.path.join(OUTPUT_DIR, "geom_info.json")
    with open(global_info_path, "w", encoding='utf-8') as f:
        json.dump(global_geom_info, f, indent=4)

if __name__ == "__main__":
    organize_dataset()
    print(f"\nGotowe! Pliki zostały uporządkowane w folderze {OUTPUT_DIR}.")