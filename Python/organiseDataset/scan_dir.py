import os
from pathlib import Path

def scan_dir(root_path):
    """Skanuje katalog i zwraca listę wszystkich plików z ich rozmiarami."""
    root = Path(root_path)
    files_info = []
    
    for path in root.rglob('*'):
        if path.is_file():
            size_mb = path.stat().st_size / (1024 * 1024)
            files_info.append((path, size_mb))
    
    return files_info


def format_data_path(path: Path) -> str:
    parts = path.parts
    for index, part in enumerate(parts):
        if part.lower() == 'data':
            return Path(*parts[index:]).as_posix()
    return path.as_posix()


def main():
    root_path = r"C:\Users\Admin\Documents\Projects\ecg project\Scripts\data\Dataset"
    files_info = scan_dir(root_path)

    with open("organiseDataset/files_report.txt", "w") as report:
        for path, size in files_info:
            report.write(f"{format_data_path(path)}: {size:.2f} MB\n")
    


if __name__ == "__main__":
    main()