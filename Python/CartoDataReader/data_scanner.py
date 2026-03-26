import os
import re
from collections import defaultdict

from pathlib import Path
from dotenv import load_dotenv
load_dotenv()
data_path = Path(os.getenv("ENV_DATA_PATH")).resolve()



def get_pattern(filename):
    """Zamienia cyfry w nazwie na '#' i usuwa rozszerzenie, by pogrupować pliki."""
    name = Path(filename).stem
    # Zamień ciągi cyfr na jeden znak #
    pattern = re.sub(r'\d+', '#', name)
    return pattern + Path(filename).suffix

def smart_scan_carto(root_path):
    root = Path(root_path)
    print(f"🕵️  Inteligentne skanowanie: {root.absolute()}\n")
    
    # Grupowanie: (pattern, level) -> [lista ścieżek]
    groups = defaultdict(list)
    
    for path in root.rglob('*'):
        if path.is_file():
            lvl = len(path.relative_to(root).parts)
            pattern = get_pattern(path.name)
            groups[(pattern, lvl, path.parent)].append(path)

    print(f"{'POZIOM':<8},{'ILOŚĆ':<6},{'ŚREDNI ROZMIAR [MB]':<15},{'WZÓR NAZWY / ŚCIEŻKA'}")
    print("-" * 90)

    # Sortujemy po poziomie i rozmiarze
    sorted_groups = sorted(groups.items(), key=lambda x: (x[0][1], -sum(f.stat().st_size for f in x[1])))

    for (pattern, lvl, folder), files in sorted_groups:
        count = len(files)
        total_size_mb = sum(f.stat().st_size for f in files) / (1024 * 1024)
        avg_size = total_size_mb / count
        
        rel_folder = folder.relative_to(root)
        
        if count > 1:
            line = f"Lvl {lvl},{count},{avg_size:.2f},{rel_folder}\\{pattern}"
        else:
            # Dla pojedynczych plików wypisz po prostu nazwę
            line = f"Lvl {lvl},{count},{total_size_mb:.2f},{rel_folder}\\{files[0].name}"
        
        # Wyświetlamy tylko jeśli plik jest "znaczący" (>0.01MB) lub jest ich dużo
        if total_size_mb > 0.1 or count > 5:
            print(line)

    print("\n" + "="*90)
    print("💡 WSKAZÓWKA: Szukaj plików z 'Points_Export.xml' (cała mapa) ")
    print("   lub ogromnych plików .txt (surowe dane).")
    print("="*90)

# smart_scan_carto('C:/Twoja/Sciezka')

# Podaj tu ścieżkę do rozpakowanych danych CARTO
smart_scan_carto(os.path.join(data_path, 'Carto3Data', 'Patient 2025_10_20', 'VT1', 'Export_VT1-02_02_2026-12-24-27'))
