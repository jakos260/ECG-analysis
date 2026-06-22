"""
Main script to organize patient ECG dataset.
Loops through patient directories and applies selected functions.
"""
import os
import re
import json
from pathlib import Path
from organising_functions import (
    copy_patient_files,
    clean_file_names,
    create_info_json,
    get_patient_data,
    create_model_xml,
    move_leads_loc_files,
    move_leads_loc_files_back,
    move_lead_files_from_model,
    remove_empty_dirs,
    strip_file_prefix,
    organise_dataset_from_IKEM
)


# ============================================================================
# CONFIGURATION
# ============================================================================

ROOT_DIR = r"C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data"
BASE_DIR = os.path.join(ROOT_DIR, "raw")
MODELS_DIR = os.path.join(BASE_DIR, "Models")
MAPPER_DIR = os.path.join(BASE_DIR, "Mapper")
ECG_DATA_DIR = os.path.join(MAPPER_DIR, "ECG_DATA")
OUTPUT_DIR = os.path.join(ROOT_DIR, "Dataset")

# List of processing steps to execute (order matters!)
# Add/remove step names as needed - just update the list
PROCESSING_STEPS = [    
    # "clean_file_names",
    # "remove_empty_dirs"
    # "move_lead_files_from_model",
    # "move_leads_loc_files_back",
    # "move_leads_loc_files",
    # "copy_files",
    # "create_info_json",
    # "create_model_xml",
    # "strip_file_prefix",
    "organise_dataset_from_IKEM"
]

# ============================================================================
# STEP REGISTRY - Maps step names to their implementation
# ============================================================================
# Add new steps here by defining the function implementation
# Format: "step_name": (step_index, description, function_call)

STEP_REGISTRY = {
    "copy_files": {
        "index": 1,
        "description": "Copying patient files",
        "func": lambda model_path, clean_name, core_id, dest_dir: copy_patient_files(
            model_path,
            dest_dir,
            MAPPER_DIR,
            ECG_DATA_DIR,
            core_id
        )
    },
    "create_info_json": {
        "index": 2,
        "description": "Creating info.json",
        "func": lambda model_path, clean_name, core_id, dest_dir: create_info_json(dest_dir)
    },
    "clean_file_names": {
        "index": 3,
        "description": "Cleaning file names",
        "func": lambda model_path, clean_name, core_id, dest_dir: clean_file_names(
            dest_dir,
            [f"{clean_name}_", f"{core_id}_"]
        )
    },
    "create_model_xml": {
        "index": 4,
        "description": "Creating model XML",
        "func": lambda model_path, clean_name, core_id, dest_dir: create_model_xml(Path(dest_dir), clean_name)
    },
    "move_leads_loc_files": {
        "index": 5,
        "description": "Moving leads location files",
        "func": lambda model_path, clean_name, core_id, dest_dir: move_leads_loc_files(dest_dir)
    },
    "move_leads_loc_files_back": {
        "index": 6,
        "description": "Moving leads location files back to ecgs",
        "func": lambda model_path, clean_name, core_id, dest_dir: move_leads_loc_files_back(dest_dir)
    },
    "move_lead_files_from_model": {
        "index": 7,
        "description": "Moving .lead and .lead12 files from model to leads",
        "func": lambda model_path, clean_name, core_id, dest_dir: move_lead_files_from_model(dest_dir)
    },
    "remove_empty_dirs": {
        "index": 8,
        "description": "Removing empty directories from patient folder",
        "func": lambda model_path, clean_name, core_id, dest_dir: remove_empty_dirs(dest_dir)
    },
    "strip_file_prefix": {
        "index": 9,
        "description": "Stripping file prefixes",
        "func": lambda model_path, clean_name, core_id, dest_dir: strip_file_prefix(Path(dest_dir) / "model", "model")
    },
    "organise_dataset_from_IKEM": {
        "index": 10,
        "description": "Organising dataset from IKEM",
        "func": lambda model_path, clean_name, core_id, dest_dir: organise_dataset_from_IKEM(ROOT_DIR, BASE_DIR, MODELS_DIR, MAPPER_DIR, ECG_DATA_DIR, OUTPUT_DIR)
    }
}


def extract_patient_id(patient_folder: str) -> str:
    """Extract patient core ID from folder name."""
    core_id_match = re.search(r'(Pat\d{3})', patient_folder)
    return core_id_match.group(1) if core_id_match else None


def get_clean_patient_name(patient_folder: str) -> str:
    """Get clean patient name by removing '_model' suffix if present."""
    clean_name = patient_folder
    if clean_name.endswith("_model"):
        clean_name = clean_name[:-6]
    return clean_name


def process_patient(patient_folder: str, model_path: str, output_dir: str) -> None:
    """
    Process a single patient directory.
    
    Args:
        patient_folder: Name of patient folder
        model_path: Full path to patient model directory
        output_dir: Output directory for organized data
    """
    clean_patient_name = get_clean_patient_name(patient_folder)
    core_id = extract_patient_id(patient_folder)

    print(f"\n{'='*70}")
    print(f"Processing: {patient_folder} -> {clean_patient_name}")
    print(f"{'='*70}")

    dest_patient_dir = os.path.join(output_dir, clean_patient_name)

    # Execute registered processing steps
    for step_name in PROCESSING_STEPS:
        if step_name not in STEP_REGISTRY:
            print(f"  ✗ Unknown step: {step_name}")
            continue

        step_info = STEP_REGISTRY[step_name]
        index = step_info["index"]
        description = step_info["description"]
        func = step_info["func"]

        print(f"  [{index}/{len(PROCESSING_STEPS)}] {description}...")
        try:
            func(model_path, clean_patient_name, core_id, dest_patient_dir)
            print(f"       ✓ Success")
        except Exception as e:
            print(f"       ✗ Error: {e}")


def collect_patient_data(output_dir: str) -> None:
    """
    Collect all patient data and save to geom_info.json.
    
    Args:
        output_dir: Output directory containing patient folders
    """
    print("\n" + "="*70)
    print("Collecting patient data into geom_info.json...")
    print("="*70)

    global_geom_info = {"patients": []}

    for patient_folder in sorted(os.listdir(output_dir)):
        patient_path = os.path.join(output_dir, patient_folder)
        if not os.path.isdir(patient_path):
            continue

        try:
            patient_data = get_patient_data(patient_path, patient_folder)
            global_geom_info["patients"].append(patient_data)
            print(f"  ✓ {patient_folder}")
        except Exception as e:
            print(f"  ✗ {patient_folder}: {e}")

    # Save global info
    global_info_path = os.path.join(output_dir, "geom_info.json")
    with open(global_info_path, "w", encoding='utf-8') as f:
        json.dump(global_geom_info, f, indent=4)
    
    print(f"\n✓ Saved: {global_info_path}")


def main():
    """Main execution function."""
    destination_only_steps = {"move_leads_loc_files", "move_leads_loc_files_back", "move_lead_files_from_model", "remove_empty_dirs", "strip_file_prefix"}
    dataset_wide_steps = {"organise_dataset_from_IKEM"}
    dataset_only_mode = len(PROCESSING_STEPS) > 0 and all(step in destination_only_steps for step in PROCESSING_STEPS)
    dataset_wide_mode = len(PROCESSING_STEPS) > 0 and all(step in dataset_wide_steps for step in PROCESSING_STEPS)

    if dataset_only_mode:
        if not os.path.exists(OUTPUT_DIR):
            print(f"ERROR: Dataset directory not found: {OUTPUT_DIR}")
            return
    else:
        # Create output directory if not exists
        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)

        # Check if source directory exists
        if not os.path.exists(MODELS_DIR):
            print(f"ERROR: Models directory not found: {MODELS_DIR}")
            return

    print("\n" + "="*70)
    print("DATASET ORGANIZATION SCRIPT")
    print("="*70)
    print(f"Models directory: {MODELS_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"\nActive steps ({len(PROCESSING_STEPS)}):")
    for i, step in enumerate(PROCESSING_STEPS, 1):
        if step in STEP_REGISTRY:
            desc = STEP_REGISTRY[step]["description"]
            print(f"  [{i}] {step}: {desc}")
        else:
            print(f"  [{i}] {step}: [UNKNOWN]")

    # Process each patient
    patient_count = 0
    if dataset_wide_mode:
        step_info = STEP_REGISTRY[PROCESSING_STEPS[0]]
        print(f"\nRunning dataset-wide step: {PROCESSING_STEPS[0]}...")
        try:
            step_info["func"](None, None, None, OUTPUT_DIR)
            print("       ✓ Success")
        except Exception as e:
            print(f"       ✗ Error: {e}")
        patient_count = sum(
            1 for folder_name in os.listdir(MODELS_DIR)
            if os.path.isdir(os.path.join(MODELS_DIR, folder_name)) and folder_name.startswith("IKEM_")
        )
    elif dataset_only_mode:
        for patient_folder in sorted(os.listdir(OUTPUT_DIR)):
            patient_path = os.path.join(OUTPUT_DIR, patient_folder)
            if not os.path.isdir(patient_path):
                continue

            process_patient(patient_folder, patient_path, OUTPUT_DIR)
            patient_count += 1
    else:
        for patient_folder in sorted(os.listdir(MODELS_DIR)):
            model_path = os.path.join(MODELS_DIR, patient_folder)
            if not os.path.isdir(model_path):
                continue

            process_patient(patient_folder, model_path, OUTPUT_DIR)
            patient_count += 1

    # Collect all patient data
    # if patient_count > 0:
    #     collect_patient_data(OUTPUT_DIR)

    print("\n" + "="*70)
    print(f"✓ COMPLETE! Processed {patient_count} patients.")
    print(f"  Output: {OUTPUT_DIR}")
    print("="*70)


if __name__ == "__main__":
    main()
