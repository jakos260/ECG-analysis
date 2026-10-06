"""Plan and apply reviewed Dataset reorganization instructions."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil
import stat


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_ROOT = SCRIPT_DIR.parents[1] / "data"
DEFAULT_INSTRUCTION_FILE = SCRIPT_DIR / "dataset_rebuild_instructions.json"
PATIENT_PATTERN = re.compile(r"Pat\d{3}", re.IGNORECASE)
TIMESTAMP_PATTERN = re.compile(
    r"\d{4}-\d{1,2}-\d{1,2}-\d{1,2}-\d{1,2}-\d{1,2}"
)
BASELINE_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:baseline|baeline|bsl)(?:[-_ ]?([12])(?=$|[^0-9]))?",
    re.IGNORECASE,
)
MEDIAN_BSM_PATTERN = re.compile(r".*\.bsm_\d\.[a-z]medianecg", re.IGNORECASE)
MODEL_WORD_PATTERN = re.compile(r"model", re.IGNORECASE)
SIGNAL_FAMILY_SUFFIX_PATTERN = re.compile(
    r"\.(?:bsm|ecg)(?:_\d)?(?:\.(?:medianecg|amedianecg|pmedianecg|vmedianecg))?$",
    re.IGNORECASE,
)


def normalize_baseline_name(name: str) -> str:
    """Normalize Baseline/BSL number tokens while preserving the rest of a name."""
    return BASELINE_PATTERN.sub(
        lambda match: f"BSL{match.group(1) or '1'}", name
    )


def normalize_ventricle_name(name: str) -> str:
    """Pluralize singular ventricle tokens without changing ventricles already plural."""
    return re.sub(
        r"ventricle(?!s)",
        lambda match: match.group() + "s",
        name,
        flags=re.IGNORECASE,
    )


def normalize_destination_name(name: str) -> str:
    return normalize_ventricle_name(name.replace("(", "").replace(")", ""))


def patient_id_from_name(name: str) -> str | None:
    match = PATIENT_PATTERN.search(name)
    return match.group().upper() if match else None


def recording_timestamp(name: str) -> str | None:
    match = TIMESTAMP_PATTERN.search(name)
    return match.group() if match else None


def normalize_timestamp(timestamp: str) -> str:
    year, month, day, hour, minute, second = (int(part) for part in timestamp.split("-"))
    return f"{year:04d}-{month:02d}-{day:02d}-{hour:02d}-{minute:02d}-{second:02d}"


def recording_key(name: str, is_source: bool = False) -> str | None:
    stem = Path(name).stem
    stem = TIMESTAMP_PATTERN.sub("", stem).strip(" _-")
    if is_source:
        stem = re.sub(r"^Pat\d{3}_nr_", "", stem, flags=re.IGNORECASE)

    baseline_match = BASELINE_PATTERN.search(stem)
    if baseline_match:
        return f"BSL{baseline_match.group(1) or '1'}"

    if is_source:
        configuration_match = re.match(r"Conf[-_ ]?(\d+)", stem, re.IGNORECASE)
        if configuration_match:
            return configuration_match.group(1)

    number_match = re.match(r"(\d+)", stem)
    return number_match.group(1) if number_match else None


def compact_recording_stem(name: str, is_source: bool = False) -> str | None:
    timestamp = recording_timestamp(name)
    key = recording_key(name, is_source=is_source)
    if not timestamp or not key:
        return None
    return f"{key}_{normalize_timestamp(timestamp)}"


def canonical_recording_label(name: str, is_source: bool = False) -> str:
    stem = Path(name).stem
    stem = TIMESTAMP_PATTERN.sub("", stem)
    if is_source:
        stem = re.sub(r"^Pat\d{3}_nr_", "", stem, flags=re.IGNORECASE)
        stem = re.sub(r"^Conf[-_ ]?(\d+)[ _-]*", r"\1_", stem, flags=re.IGNORECASE)
    stem = normalize_baseline_name(stem)
    return re.sub(r"[^a-z0-9]+", "", stem.lower())


def source_fallback_stem(source_name: str) -> str:
    stem = compact_recording_stem(source_name, is_source=True)
    return stem or normalize_baseline_name(Path(source_name).stem)


def bsm_index(subject_root: Path, patient_id: str) -> dict[tuple[str, str], list[Path]]:
    index: dict[tuple[str, str], list[Path]] = {}
    data_dir = subject_root / "BSM" / "ECG_DATA"
    if not data_dir.is_dir():
        return index

    for path in data_dir.iterdir():
        if not path.is_file() or path.suffix.lower() != ".bsm":
            continue
        timestamp = recording_timestamp(path.name)
        if timestamp:
            index.setdefault((patient_id, normalize_timestamp(timestamp)), []).append(path)
    return index


def match_bsm(source_path: Path, index: dict[tuple[str, str], list[Path]]) -> tuple[Path | None, str]:
    patient_id = patient_id_from_name(source_path.name)
    timestamp = recording_timestamp(source_path.name)
    if not patient_id or not timestamp:
        return None, "source filename has no patient ID or recording timestamp"

    candidates = index.get((patient_id, normalize_timestamp(timestamp)), [])
    if len(candidates) == 1:
        return candidates[0], "matched by patient ID and timestamp"
    if not candidates:
        return None, "no BSM has the same patient ID and timestamp"

    source_label = canonical_recording_label(source_path.name, is_source=True)
    matching = [
        candidate
        for candidate in candidates
        if canonical_recording_label(candidate.name) == source_label
    ]
    if len(matching) == 1:
        return matching[0], "timestamp collision resolved by normalized recording label"
    return None, "multiple BSM files share the timestamp and labels do not resolve the match"


def operation(op: str, **values: str) -> dict[str, str]:
    return {"op": op, **values}


def generate_instructions(data_root: Path) -> dict:
    data_root = data_root.resolve()
    dataset_root = data_root / "Dataset"
    source_root = data_root / "raw" / "Mapper" / "ECG_DATA"
    if not dataset_root.is_dir():
        raise FileNotFoundError(f"Dataset directory not found: {dataset_root}")
    if not source_root.is_dir():
        raise FileNotFoundError(f"Source ECG_DATA directory not found: {source_root}")

    operations: list[dict] = []
    warnings: list[str] = []
    signal_names: dict[str, dict[str, str]] = {}
    subject_by_patient: dict[str, Path] = {}

    for subject in sorted(path for path in dataset_root.iterdir() if path.is_dir()):
        patient_id = patient_id_from_name(subject.name)
        if not patient_id:
            continue
        if patient_id in subject_by_patient:
            warnings.append(
                f"Multiple subject folders contain {patient_id}: "
                f"{subject_by_patient[patient_id]} and {subject}"
            )
            continue
        subject_by_patient[patient_id] = subject

    active_mapper_roots: dict[str, Path] = {}
    bsm_indexes: dict[str, dict[tuple[str, str], list[Path]]] = {}
    for patient_id, subject in subject_by_patient.items():
        map_root = subject / "map"
        mapper_root = subject / "mapper"
        if map_root.exists() and mapper_root.exists():
            warnings.append(
                f"{subject.name}: both map and mapper exist; map was left untouched"
            )
            active_root = mapper_root
        elif map_root.exists():
            operations.append(
                operation("mv", source=str(map_root), destination=str(mapper_root))
            )
            active_root = mapper_root
        else:
            active_root = mapper_root

        active_mapper_roots[patient_id] = active_root
        index_root = map_root if map_root.is_dir() else mapper_root
        bsm_indexes[patient_id] = bsm_index(index_root, patient_id)

        current_root = map_root if map_root.is_dir() else mapper_root
        if not current_root.is_dir():
            continue

        for current_path in sorted(path for path in current_root.rglob("*") if path.is_file()):
            relative_path = current_path.relative_to(current_root)
            target_path = active_root / relative_path
            lower_name = current_path.name.lower()

            if lower_name.endswith(".iecg"):
                operations.append(operation("rm", destination=str(target_path)))
                continue

            if len(relative_path.parts) == 1 and lower_name.endswith((".imap", ".imaplog")):
                normalized = "".join(character for character in lower_name if character.isalnum())
                group = "12ECG" if "12ecg" in normalized or "ecg12" in normalized else "BSM"
                new_name = normalize_baseline_name(current_path.name)
                destination = active_root / group / new_name
                operations.append(
                    operation("mv", source=str(target_path), destination=str(destination))
                )
                continue

            if len(relative_path.parts) == 1 and (
                lower_name.endswith((".bsm", ".ecg")) or MEDIAN_BSM_PATTERN.fullmatch(lower_name)
            ):
                group = "12ECG" if lower_name.endswith(".ecg") else "BSM"
                data_destination = active_root / group / "ECG_DATA" / normalize_baseline_name(current_path.name)
                operations.append(
                    operation("mv", source=str(target_path), destination=str(data_destination))
                )
                continue

            normalized_name = normalize_baseline_name(current_path.name)
            if normalized_name != current_path.name:
                destination = active_root / relative_path.parent / normalized_name
                operations.append(
                    operation("mv", source=str(target_path), destination=str(destination))
                )

    source_ecg_files = sorted(
        path for path in source_root.iterdir()
        if path.is_file() and path.suffix.lower() == ".ecg"
    )
    names_by_subject: dict[str, dict[str, str]] = {}
    reserved_copy_destinations: set[Path] = set()

    for source_path in source_ecg_files:
        patient_id = patient_id_from_name(source_path.name)
        subject = subject_by_patient.get(patient_id or "")
        if not patient_id or not subject:
            warnings.append(f"No Dataset subject found for source file: {source_path.name}")
            continue

        bsm_match, match_reason = match_bsm(source_path, bsm_indexes[patient_id])
        if bsm_match:
            destination_stem = compact_recording_stem(bsm_match.name)
            if not destination_stem:
                destination_stem = source_fallback_stem(source_path.name)
                warnings.append(
                    f"{source_path.name}: matched BSM name has no recording key or timestamp; "
                    f"using {destination_stem}.ecg"
                )
        else:
            destination_stem = source_fallback_stem(source_path.name)
            warnings.append(f"{source_path.name}: {match_reason}; using {destination_stem}.ecg")

        destination_name = destination_stem + ".ecg"
        mapper_root = active_mapper_roots[patient_id]
        destination = mapper_root / "12ECG" / "ECG_DATA" / destination_name
        if destination in reserved_copy_destinations:
            warnings.append(f"Copy skipped because another source maps to: {destination}")
            continue

        prior_name = names_by_subject.setdefault(subject.name, {})
        destination_key = Path(destination_name).stem
        if destination_key in prior_name:
            warnings.append(
                f"Copy skipped because {subject.name} maps multiple sources to {destination_key}"
            )
            continue

        reserved_copy_destinations.add(destination)
        prior_name[destination_key] = source_path.stem
        if destination.exists():
            warnings.append(f"Copy skipped because destination already exists: {destination}")
        else:
            operations.append(
                operation("cp", source=str(source_path), destination=str(destination))
            )

    signal_names = names_by_subject
    for subject_name, names in sorted(signal_names.items()):
        subject_root = dataset_root / subject_name
        mapper_root = active_mapper_roots[patient_id_from_name(subject_name) or ""]
        operations.append(
            {
                "op": "write_json",
                "path": str(mapper_root / "signal_names.json"),
                "data": names,
                "overwrite": True,
            }
        )

    operation_counts = {
        op: sum(item["op"] == op for item in operations)
        for op in ("mv", "cp", "rm", "write_json")
    }
    for item in operations:
        for key in ("source", "destination", "path"):
            if key in item:
                item[key] = str(Path(item[key]).relative_to(data_root))

    return {
        "name": "Dataset Cartographer 2.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "data_root": str(data_root),
        "operations": operations,
        "signal_names": signal_names,
        "warnings": warnings,
        "summary": {
            "subjects_scanned": len(subject_by_patient),
            "source_ecg_files": len(source_ecg_files),
            "operations": operation_counts,
        },
    }


def cleaned_model_name(file_name: str, source_subject_name: str, patient_id: str | None) -> str:
    subject_prefix = re.sub(r"_model$", "", source_subject_name, flags=re.IGNORECASE)
    subject_parts = subject_prefix.split("_")
    prefix_candidates = {source_subject_name, subject_prefix}
    prefix_candidates.update(
        "_".join(subject_parts[index:])
        for index in range(1, len(subject_parts))
    )
    cleaned = file_name
    for prefix in sorted(prefix_candidates, key=len, reverse=True):
        cleaned_candidate = re.sub(
            rf"^{re.escape(prefix)}(?=$|[ _.-])",
            "",
            cleaned,
            count=1,
            flags=re.IGNORECASE,
        )
        if cleaned_candidate != cleaned:
            cleaned = cleaned_candidate
            break
    if patient_id:
        cleaned = re.sub(re.escape(patient_id), "", cleaned, flags=re.IGNORECASE)
    cleaned = MODEL_WORD_PATTERN.sub("", cleaned)

    if cleaned.casefold() == ".xml":
        cleaned = "subject" + cleaned
    else:
        cleaned = cleaned.lstrip(" ._-").rstrip(" _-")
    return normalize_destination_name(cleaned or "subject_file")


def mapper_signal_destination(file_name: str) -> tuple[str, str] | None:
    patient_id = patient_id_from_name(file_name)
    timestamp_match = TIMESTAMP_PATTERN.search(file_name)
    key = recording_key(file_name, is_source=True)
    if not patient_id or not timestamp_match or not key:
        return None

    suffix = file_name[timestamp_match.end():]
    if suffix.lower().startswith(".bsm"):
        group = "BSM"
    elif suffix.lower().startswith(".ecg"):
        group = "12ECG"
    else:
        return None

    return group, normalize_destination_name(key + suffix)


def metadata_group(file_name: str) -> str:
    normalized = "".join(character for character in file_name.lower() if character.isalnum())
    return "12ECG" if "12ecg" in normalized or "ecg12" in normalized else "BSM"


def signal_name_without_extension(file_name: str) -> str:
    return SIGNAL_FAMILY_SUFFIX_PATTERN.sub("", Path(file_name).name)


def cleaned_metadata_name(file_name: str, patient_id: str) -> str:
    cleaned = re.sub(re.escape(patient_id), "", file_name, flags=re.IGNORECASE)
    return normalize_destination_name(normalize_baseline_name(cleaned.lstrip(" _-")))


def generate_rebuild_instructions(data_root: Path) -> dict:
    data_root = data_root.resolve()
    models_root = data_root / "raw" / "Models"
    mapper_root = data_root / "raw" / "Mapper"
    dataset_root = data_root / "Dataset"
    if not models_root.is_dir():
        raise FileNotFoundError(f"Models directory not found: {models_root}")
    if not mapper_root.is_dir():
        raise FileNotFoundError(f"Mapper directory not found: {mapper_root}")

    model_subjects = sorted(
        (path for path in models_root.iterdir() if path.is_dir()),
        key=lambda path: path.name.casefold(),
    )
    if not model_subjects:
        raise ValueError(f"No subject model directories found: {models_root}")

    operations: list[dict] = [operation("rm", destination=str(dataset_root))]
    warnings: list[str] = []
    subject_registry: dict[str, dict[str, str]] = {}
    subjects_by_patient: dict[str, tuple[str, Path, str]] = {}
    reserved_destinations: set[Path] = set()
    signal_names: dict[str, dict[str, dict[str, str]]] = {}

    for index, source_subject in enumerate(model_subjects, start=1):
        subject_name = f"subject_{index:03d}"
        original_name = source_subject.name
        display_name = normalize_destination_name(
            re.sub(r"_model$", "", original_name, flags=re.IGNORECASE)
        )
        patient_id = patient_id_from_name(original_name)
        subject_registry[subject_name] = {
            "subject_name": display_name,
            "source_model_dir": str(source_subject.relative_to(data_root)),
        }
        if patient_id:
            if patient_id in subjects_by_patient:
                warnings.append(
                    f"Duplicate model directories for {patient_id}: "
                    f"{subjects_by_patient[patient_id][1]} and {source_subject}"
                )
            else:
                subjects_by_patient[patient_id] = (subject_name, source_subject, original_name)

        operations.append(
            {
                "op": "write_json",
                "path": str(dataset_root / subject_name / "subject.json"),
                "data": {"subject_name": display_name},
            }
        )

        model_files = sorted(path for path in source_subject.rglob("*") if path.is_file())
        for source_path in model_files:
            relative_parts = list(source_path.relative_to(source_subject).parts)
            if relative_parts[0].casefold() == "ecgs":
                relative_parts.pop(0)
                cleaned_parts = [
                    cleaned_model_name(part, original_name, patient_id)
                    for part in relative_parts[:-1]
                ]
                cleaned_file = cleaned_model_name(relative_parts[-1], original_name, patient_id)
                destination = (
                    dataset_root / subject_name / "signals" / "ecgs"
                    / Path(*cleaned_parts) / cleaned_file
                )
                if destination in reserved_destinations:
                    warnings.append(f"Model ecg copy skipped due to duplicate destination: {destination}")
                    continue
                reserved_destinations.add(destination)
                operations.append(
                    operation("cp", source=str(source_path), destination=str(destination))
                )
                continue

            if relative_parts and relative_parts[0].casefold() == "model":
                relative_parts.pop(0)
            cleaned_parts = [
                cleaned_model_name(part, original_name, patient_id)
                for part in relative_parts[:-1]
            ]
            cleaned_file = cleaned_model_name(relative_parts[-1], original_name, patient_id)
            destination = dataset_root / subject_name / "model" / Path(*cleaned_parts) / cleaned_file
            if destination in reserved_destinations:
                warnings.append(f"Model copy skipped due to duplicate destination: {destination}")
                continue
            reserved_destinations.add(destination)
            operations.append(
                operation("cp", source=str(source_path), destination=str(destination))
            )

    loose_model_files = sorted(path for path in models_root.iterdir() if path.is_file())
    for source_path in loose_model_files:
        patient_id = patient_id_from_name(source_path.name)
        subject_record = subjects_by_patient.get(patient_id or "")
        if not subject_record:
            warnings.append(f"Unassigned Models root file skipped: {source_path.name}")
            continue
        subject_name, _, original_name = subject_record
        destination_name = cleaned_model_name(source_path.name, original_name, patient_id)
        destination = dataset_root / subject_name / "model" / destination_name
        if destination in reserved_destinations:
            warnings.append(f"Model copy skipped due to duplicate destination: {destination}")
            continue
        reserved_destinations.add(destination)
        operations.append(
            operation("cp", source=str(source_path), destination=str(destination))
        )

    mapper_files = sorted(path for path in mapper_root.rglob("*") if path.is_file())
    mappings_by_subject: dict[str, dict[str, dict[str, str]]] = {}
    for source_path in mapper_files:
        relative_path = source_path.relative_to(mapper_root)
        patient_id = patient_id_from_name(source_path.name)
        subject_record = subjects_by_patient.get(patient_id or "")
        if not patient_id or not subject_record:
            warnings.append(f"Unassigned Mapper file skipped: {relative_path.as_posix()}")
            continue

        subject_name, _, _ = subject_record
        subject_mapper_root = dataset_root / subject_name / "mapper"
        subject_signals_root = dataset_root / subject_name / "signals"
        if source_path.suffix.lower() == ".iecg":
            destination_name = cleaned_metadata_name(source_path.name, patient_id)
            destination = subject_signals_root / destination_name
        elif relative_path.parts and relative_path.parts[0].casefold() == "ecg_data":
            classified = mapper_signal_destination(source_path.name)
            if not classified:
                warnings.append(f"Unrecognized Mapper ECG_DATA file skipped: {relative_path.as_posix()}")
                continue

            group, destination_name = classified
            destination = subject_signals_root / "ECG_DATA" / group / destination_name
            mapping_group = mappings_by_subject.setdefault(
                subject_name, {"BSM": {}, "12ECG": {}}
            )[group]
            new_name = signal_name_without_extension(destination_name)
            old_name = signal_name_without_extension(source_path.name)
            if new_name in mapping_group and mapping_group[new_name] != old_name:
                warnings.append(
                    f"signal_names collision for {subject_name}/{group}/{new_name}; "
                    f"file skipped: {relative_path.as_posix()}"
                )
                continue
            mapping_group[new_name] = old_name
        else:
            destination_name = cleaned_metadata_name(source_path.name, patient_id)
            destination = subject_mapper_root / destination_name

        if destination in reserved_destinations:
            warnings.append(f"Mapper copy skipped due to duplicate destination: {destination}")
            continue
        reserved_destinations.add(destination)
        operations.append(
            operation("cp", source=str(source_path), destination=str(destination))
        )

    for subject_name, names_by_group in sorted(mappings_by_subject.items()):
        operations.append(
            {
                "op": "write_json",
                "path": str(dataset_root / subject_name / "signals" / "signal_names.json"),
                "data": names_by_group,
            }
        )

    operation_counts = {
        op: sum(item["op"] == op for item in operations)
        for op in ("rm", "cp", "write_json")
    }
    for item in operations:
        for key in ("source", "destination", "path"):
            if key in item:
                item[key] = str(Path(item[key]).relative_to(data_root))

    model_file_count = sum(
        1 for subject in model_subjects for path in subject.rglob("*") if path.is_file()
    ) + len(loose_model_files)
    return {
        "name": "Dataset Rebuild Instructions",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "data_root": str(data_root),
        "subject_registry": subject_registry,
        "signal_names": mappings_by_subject,
        "operations": operations,
        "warnings": warnings,
        "summary": {
            "subjects": len(model_subjects),
            "model_files_scanned": model_file_count,
            "mapper_files_scanned": len(mapper_files),
            "operations": operation_counts,
        },
    }


def write_instruction_file(data_root: Path, output_path: Path) -> dict:
    instructions = generate_rebuild_instructions(data_root)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(instructions, indent=2), encoding="utf-8")
    return instructions


def ensure_within_root(path: Path, data_root: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(data_root.resolve()):
        raise ValueError(f"Instruction path escapes data root: {path}")
    return resolved


def resolve_instruction_path(path_value: str, data_root: Path) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        raise ValueError(f"Instruction paths must be relative to data_root: {path_value}")
    return ensure_within_root(data_root / path, data_root)


def remove_file(path: Path) -> None:
    try:
        path.unlink()
    except PermissionError:
        path.chmod(stat.S_IWRITE | stat.S_IREAD)
        path.unlink()


def apply_instructions(instruction_path: Path, dry_run: bool = False) -> int:
    instructions = json.loads(instruction_path.read_text(encoding="utf-8"))
    data_root = Path(instructions["data_root"]).resolve()
    operations = instructions["operations"]

    for item in operations:
        op = item["op"]
        if op in ("mv", "cp"):
            source = resolve_instruction_path(item["source"], data_root)
            destination = resolve_instruction_path(item["destination"], data_root)
            print(f"{op} {source} -> {destination}")
            if dry_run:
                continue
            if not source.exists():
                raise FileNotFoundError(f"Instruction source not found: {source}")
            if destination.exists():
                raise FileExistsError(f"Instruction destination already exists: {destination}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            if op == "mv":
                shutil.move(str(source), str(destination))
            else:
                if not source.is_file():
                    raise ValueError(f"Copy source is not a file: {source}")
                shutil.copy2(source, destination)
        elif op == "rm":
            destination = resolve_instruction_path(item["destination"], data_root)
            print(f"rm {destination}")
            if dry_run:
                continue
            if destination == (data_root / "Dataset").resolve():
                if destination.exists():
                    if not destination.is_dir():
                        raise NotADirectoryError(f"Dataset removal target is not a directory: {destination}")

                    def make_writable_and_retry(func, path, exc_info):
                        os.chmod(path, stat.S_IWRITE | stat.S_IREAD | stat.S_IXUSR)
                        func(path)

                    shutil.rmtree(destination, onerror=make_writable_and_retry)
            elif destination.is_file():
                remove_file(destination)
            else:
                raise FileNotFoundError(f"Removal target is not a file: {destination}")
        elif op == "write_json":
            destination = resolve_instruction_path(item["path"], data_root)
            print(f"write_json {destination}")
            if dry_run:
                continue
            if destination.exists() and not item.get("overwrite", False):
                raise FileExistsError(f"JSON destination already exists: {destination}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(json.dumps(item["data"], indent=2), encoding="utf-8")
        else:
            raise ValueError(f"Unknown instruction operation: {op}")

    return len(operations)


def main() -> None:
    parser = argparse.ArgumentParser(description="Dataset Cartographer 2.0")
    commands = parser.add_subparsers(dest="command", required=True)

    plan_parser = commands.add_parser("plan", help="Scan files and generate instructions")
    plan_parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    plan_parser.add_argument("--output", type=Path, default=DEFAULT_INSTRUCTION_FILE)

    apply_parser = commands.add_parser("apply", help="Apply a reviewed instruction file")
    apply_parser.add_argument("instruction_file", type=Path)
    apply_parser.add_argument("--dry-run", action="store_true")

    args = parser.parse_args()
    if args.command == "plan":
        instructions = write_instruction_file(args.data_root, args.output)
        print(f"Instructions: {args.output.resolve()}")
        print(json.dumps(instructions["summary"], indent=2))
        print(f"Warnings requiring review: {len(instructions['warnings'])}")
    else:
        count = apply_instructions(args.instruction_file, dry_run=args.dry_run)
        mode = "Validated" if args.dry_run else "Applied"
        print(f"{mode} {count} instructions")


if __name__ == "__main__":
    main()