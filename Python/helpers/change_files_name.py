import os

def usun_prefix(folder_path):
    # Pobieramy samą nazwę folderu (np. "X")
    folder_name = os.path.basename(os.path.normpath(folder_path))
    
    # Przechodzimy po wszystkich plikach w folderze
    for filename in os.listdir(folder_path):
        old_path = os.path.join(folder_path, filename)
        
        # Pomijamy podfoldery
        if os.path.isdir(old_path):
            continue
        
        # Sprawdzamy, czy plik zaczyna się od prefixu "X_"
        prefix = folder_name + "_"
        if filename.startswith(prefix):
            new_filename = filename[len(prefix):]  # usuwamy prefix
            new_path = os.path.join(folder_path, new_filename)
            
            # Zmieniamy nazwę pliku
            os.rename(old_path, new_path)
            print(f"Zmieniono nazwę: {filename} → {new_filename}")

if __name__ == "__main__":
    # folder = input("Podaj ścieżkę do folderu: ").strip()
    # usun_prefix(folder)
    usun_prefix("PPD1_ECGSIM_NEW")