"""
Utility functions for organizing patient ECG dataset.
"""
import os
import shutil
import json
from pathlib import Path
from typing import Dict, Tuple, Optional


def _normalize_ikem_patient_name(folder_name: str) -> str:
    """Convert an IKEM model folder name into the dataset patient folder name."""
    patient_name = folder_name
    if patient_name.endswith("_model"):
        patient_name = patient_name[:-len("_model")]
    return patient_name


def _strip_prefix_from_filename(file_name: str, prefixes: list) -> str:
    """Remove the first matching prefix from a filename and trim separator characters."""
    for prefix in sorted(prefixes, key=len, reverse=True):
        if file_name.startswith(prefix):
            return file_name[len(prefix):].lstrip(" _-")
    return file_name


def _normalize_copied_filename(file_name: str) -> str:
    """Normalize a copied filename for the output dataset."""
    return file_name.replace(" ", "_")


def clean_patient_data(dest_patient_dir: str) -> None:
    """Remove an existing patient dataset directory before rebuilding it."""
    if os.path.exists(dest_patient_dir):
        shutil.rmtree(dest_patient_dir)
    os.makedirs(dest_patient_dir, exist_ok=True)


def copy_model_files(src_dir: str, dest_patient_dir: str, patient_name: str) -> None:
    """Copy IKEM model files into model/leads folders with prefixes removed."""
    dest_model_dir = os.path.join(dest_patient_dir, "model")
    dest_leads_dir = os.path.join(dest_patient_dir, "leads")
    os.makedirs(dest_model_dir, exist_ok=True)
    os.makedirs(dest_leads_dir, exist_ok=True)

    prefix_candidates = [
        f"{patient_name}_model_",
        f"{patient_name}_model",
        f"IKEM_{patient_name}_model",
        f"IKEM_{patient_name}_model_",
    ]

    for item in os.listdir(src_dir):
        src_item = os.path.join(src_dir, item)
        if not os.path.isfile(src_item):
            continue

        cleaned_name = _strip_prefix_from_filename(item, prefix_candidates)
        if cleaned_name == ".sternummarker":
            cleaned_name = "model.sternummarker"
        lower_name = cleaned_name.lower()

        if lower_name.endswith(".lead") or lower_name.endswith(".lead12"):
            shutil.copy2(src_item, os.path.join(dest_leads_dir, cleaned_name))
        else:
            shutil.copy2(src_item, os.path.join(dest_model_dir, cleaned_name))


def copy_mapper_files(mapper_dir: str, ecg_data_dir: str, dest_patient_dir: str, patient_core_id: str) -> None:
    """Copy IKEM mapper files into signals/map folders with prefixes removed."""
    dest_signals_dir = os.path.join(dest_patient_dir, "signals")
    dest_map_dir = os.path.join(dest_patient_dir, "map")
    os.makedirs(dest_signals_dir, exist_ok=True)
    os.makedirs(dest_map_dir, exist_ok=True)

    signal_prefixes = [f"{patient_core_id}_nr_Conf-", f"{patient_core_id}_nr_"]
    map_prefixes = [f"{patient_core_id}_"]
    signal_suffixes = (".bsm.medianecg", ".bsm", ".bms", ".ecg")
    map_suffixes = (".iecg", ".imap", ".imaplog")

    if os.path.exists(ecg_data_dir):
        for file_name in os.listdir(ecg_data_dir):
            src_path = os.path.join(ecg_data_dir, file_name)
            if not os.path.isfile(src_path):
                continue
            if not any(file_name.startswith(prefix) for prefix in signal_prefixes):
                continue

            lower_name = file_name.lower()
            if lower_name.endswith(signal_suffixes):
                cleaned_name = _strip_prefix_from_filename(file_name, signal_prefixes)
                cleaned_name = _normalize_copied_filename(cleaned_name)
                shutil.copy2(src_path, os.path.join(dest_signals_dir, cleaned_name))

    if os.path.exists(mapper_dir):
        for file_name in os.listdir(mapper_dir):
            src_path = os.path.join(mapper_dir, file_name)
            if not os.path.isfile(src_path):
                continue
            if not any(file_name.startswith(prefix) for prefix in map_prefixes):
                continue

            lower_name = file_name.lower()
            if lower_name.endswith(map_suffixes):
                cleaned_name = _strip_prefix_from_filename(file_name, map_prefixes)
                cleaned_name = _normalize_copied_filename(cleaned_name)
                shutil.copy2(src_path, os.path.join(dest_map_dir, cleaned_name))


def organise_dataset_from_IKEM(root_dir: str, base_dir: str, models_dir: str, mapper_dir: str, ecg_data_dir: str, output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)

    for folder_name in sorted(os.listdir(models_dir)):
        src_model_dir = os.path.join(models_dir, folder_name)
        if not os.path.isdir(src_model_dir):
            continue
        if not folder_name.startswith("IKEM_"):
            continue

        patient_name = _normalize_ikem_patient_name(folder_name)
        patient_core_id = next((part for part in patient_name.split("_") if part.startswith("Pat")), patient_name)
        dest_patient_dir = os.path.join(output_dir, patient_name)

        clean_patient_data(dest_patient_dir)
        copy_model_files(src_model_dir, dest_patient_dir, patient_name)
        copy_mapper_files(mapper_dir, ecg_data_dir, dest_patient_dir, patient_core_id)

def copy_patient_files(src_dir: str, dest_patient_dir: str, mapper_dir: str, 
                       ecg_data_dir: str, core_id: Optional[str]) -> None:
    """
    Copy patient model files and ECG data to destination.
    
    Args:
        src_dir: Source patient directory
        dest_patient_dir: Destination patient directory
        mapper_dir: Mapper directory containing ECG files
        ecg_data_dir: ECG data directory
        core_id: Patient ID for finding related files
    """
    dest_model_dir = os.path.join(dest_patient_dir, "model")
    dest_ecgs_dir = os.path.join(dest_patient_dir, "ecgs")
    dest_ventricular_dir = os.path.join(dest_patient_dir, "ventricular_beats")

    os.makedirs(dest_model_dir, exist_ok=True)
    os.makedirs(dest_ecgs_dir, exist_ok=True)

    # Copy files preserving structure
    for item in os.listdir(src_dir):
        src_item = os.path.join(src_dir, item)
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

    # Copy ECG files from Mapper directories
    if core_id:
        if os.path.exists(mapper_dir):
            for file in os.listdir(mapper_dir):
                if core_id in file and os.path.isfile(os.path.join(mapper_dir, file)):
                    shutil.copy2(os.path.join(mapper_dir, file), dest_ecgs_dir)

        if os.path.exists(ecg_data_dir):
            for file in os.listdir(ecg_data_dir):
                if core_id in file and os.path.isfile(os.path.join(ecg_data_dir, file)):
                    file_path = os.path.join(ecg_data_dir, file)
                    shutil.copy2(file_path, dest_ecgs_dir)


def clean_file_names(dest_patient_dir: str, prefixes: list) -> None:
    """
    Remove prefixes from filenames in patient directory tree.
    
    Args:
        dest_patient_dir: Patient directory
        prefixes: List of prefixes to remove
    """
    prefixes = [p for p in prefixes if p != "_"]

    for root, dirs, files in os.walk(dest_patient_dir):
        for file in files:
            new_name = file
            for prefix in prefixes:
                if new_name.startswith(prefix):
                    new_name = new_name[len(prefix):].lstrip(" _")
                    break

            if new_name != file:
                old_path = os.path.join(root, file)
                new_path = os.path.join(root, new_name)
                if not os.path.exists(new_path):
                    os.rename(old_path, new_path)


def create_info_json(dest_patient_dir: str) -> None:
    """
    Create info.json file with default values.
    
    Args:
        dest_patient_dir: Patient directory
    """
    info_data = {"qrs_begin": -1, "t_begin": -1, "BSM_leads_number": -1}
    with open(os.path.join(dest_patient_dir, "info.json"), "w", encoding='utf-8') as f:
        json.dump(info_data, f, indent=4)


def move_leads_loc_files(dest_patient_dir: str) -> None:
    """
    Move iecg, imap, and imapLog files from the ecgs directory into a leads_loc directory.

    Args:
        dest_patient_dir: Patient directory
    """
    dest_ecgs_dir = os.path.join(dest_patient_dir, "ecgs")
    leads_loc_dir = os.path.join(dest_patient_dir, "leads_loc")

    if not os.path.exists(dest_ecgs_dir):
        return

    os.makedirs(leads_loc_dir, exist_ok=True)

    for file_name in os.listdir(dest_ecgs_dir):
        lower_name = file_name.lower()
        if lower_name.endswith(".iecg") or lower_name.endswith(".imap") or lower_name.endswith(".imaplog"):
            src_path = os.path.join(dest_ecgs_dir, file_name)
            dst_path = os.path.join(leads_loc_dir, file_name)
            if os.path.isfile(src_path):
                if os.path.exists(dst_path):
                    os.remove(dst_path)
                shutil.move(src_path, dst_path)


def move_leads_loc_files_back(dest_patient_dir: str) -> None:
    """
    Move iecg, imap, and imapLog files from the leads_loc directory back into the ecgs directory.

    Args:
        dest_patient_dir: Patient directory
    """
    leads_loc_dir = os.path.join(dest_patient_dir, "leads_loc")
    dest_ecgs_dir = os.path.join(dest_patient_dir, "ecgs")

    if not os.path.exists(leads_loc_dir):
        return

    os.makedirs(dest_ecgs_dir, exist_ok=True)

    for file_name in os.listdir(leads_loc_dir):
        lower_name = file_name.lower()
        if lower_name.endswith(".iecg") or lower_name.endswith(".imap") or lower_name.endswith(".imaplog"):
            src_path = os.path.join(leads_loc_dir, file_name)
            dst_path = os.path.join(dest_ecgs_dir, file_name)
            if os.path.isfile(src_path):
                if os.path.exists(dst_path):
                    os.remove(dst_path)
                shutil.move(src_path, dst_path)


def move_lead_files_from_model(dest_patient_dir: str) -> None:
    """
    Move .lead and .lead12 files from the model directory into a leads directory.

    Args:
        dest_patient_dir: Patient directory
    """
    model_dir = os.path.join(dest_patient_dir, "model")
    leads_dir = os.path.join(dest_patient_dir, "leads")

    if not os.path.exists(model_dir):
        return

    os.makedirs(leads_dir, exist_ok=True)

    for file_name in os.listdir(model_dir):
        lower_name = file_name.lower()
        if lower_name.endswith(".lead") or lower_name.endswith(".lead12") or lower_name.endswith(".medianecg") or lower_name.endswith(".ecg")  or lower_name.endswith(".bsm"):
            src_path = os.path.join(model_dir, file_name)
            dst_path = os.path.join(leads_dir, file_name)
            if os.path.isfile(src_path):
                if os.path.exists(dst_path):
                    os.remove(dst_path)
                shutil.move(src_path, dst_path)


def remove_empty_dirs(dest_patient_dir: str) -> None:
    """
    Remove all empty directories inside a patient destination directory.

    Args:
        dest_patient_dir: Patient directory
    """
    for root, dirs, files in os.walk(dest_patient_dir, topdown=False):
        for d in dirs:
            dir_path = os.path.join(root, d)
            if os.path.isdir(dir_path) and not os.listdir(dir_path):
                os.rmdir(dir_path)


def get_patient_data(dest_patient_dir: str, clean_patient_name: str) -> Dict:
    """
    Build patient data dictionary with model and signal information.
    
    Args:
        dest_patient_dir: Patient directory
        clean_patient_name: Clean patient name
        
    Returns:
        Dictionary with patient model and signal data
    """
    dest_model_dir = os.path.join(dest_patient_dir, "model")
    dest_ecgs_dir = os.path.join(dest_patient_dir, "ecgs")
    dest_beats_dir = os.path.join(dest_patient_dir, "ventricular_beats")

    is_ecgsim = "ECGSIM" in clean_patient_name.upper()
    is_ikem = "IKEM" in clean_patient_name.upper()

    # --- MODEL SECTION ---
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

    # --- SIGNALS SECTION ---
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

    # Build measurements section with different logic for ECGSIM vs IKEM
    measurements = []

    if is_ecgsim:
        bsm_ref = ""
        ecg_ref = ""

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

        for m_name in sorted(grouped_measurements.keys()):
            measurements.append({
                "name": m_name,
                "bsm": grouped_measurements[m_name]["bsm"],
                "bsm_medianecg": grouped_measurements[m_name]["bsm_medianecg"],
                "ecg": grouped_measurements[m_name]["ecg"]
            })

    # Ground truth
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


def create_model_xml(patient_folder: Path, patient_id: str) -> None:
    """
    Create XML file for patient model.
    
    Args:
        patient_folder: Path to patient folder
        patient_id: Patient ID
    """
    xml_template = """<!DOCTYPE ANATOMY>
<ANATOMY version="1">
  <ANATOMYDATA version="1">
    <Shapes>
      <ShapeDirectory>.</ShapeDirectory>
    </Shapes>
    <patientId>{patient_id}</patientId>
    <patientHeight></patientHeight>
    <patientChest></patientChest>
    <patientAge>0</patientAge>
    <patientNotes></patientNotes>
    <anisratio>2.5</anisratio>
  </ANATOMYDATA>
</ANATOMY>"""

    model_folder = patient_folder / "model"

    if model_folder.exists() and model_folder.is_dir():
        # Rename .tri files
        for tri_file in model_folder.glob("*.tri"):
            if not tri_file.name.startswith("model_"):
                new_name = f"model_{tri_file.name}"
                new_path = model_folder / new_name
                tri_file.rename(new_path)

        # Create XML file
        xml_filename = f"{patient_id}_model.xml"
        xml_filepath = model_folder / xml_filename
        formatted_xml = xml_template.format(patient_id=patient_id)

        with open(xml_filepath, "w", encoding="utf-8") as f:
            f.write(formatted_xml)


def strip_file_prefix(target_dir: Path, prefix: str, overwrite: bool = False) -> Tuple[int, int, int]:
    """
    Strip prefix from filenames in directory.
    
    Args:
        target_dir: Target directory
        prefix: Prefix to remove
        overwrite: Whether to overwrite existing files
        
    Returns:
        Tuple of (renamed_count, skipped_count, errors_count)
    """
    prefix_with_sep = prefix + "_"
    renamed = 0
    skipped = 0
    errors = 0

    for p in sorted(target_dir.iterdir()):
        if not p.is_file():
            continue
        name = p.name
        if not name.startswith(prefix_with_sep):
            continue

        new_name = name[len(prefix_with_sep):]
        dest = target_dir / new_name

        try:
            if dest.exists():
                if overwrite:
                    dest.unlink()
                    p.rename(dest)
                    renamed += 1
                else:
                    skipped += 1
            else:
                p.rename(dest)
                renamed += 1
        except Exception as e:
            print(f"ERROR renaming {p} -> {dest}: {e}")
            errors += 1

    return renamed, skipped, errors
