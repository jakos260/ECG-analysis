classdef DatasetLoader < handle
    properties
        FolderPath      % Ścieżka do głównego folderu z podkatalogami pacjentów
        Patients        % Tablica obiektów PatientData
    end
    
    methods
        function obj = DatasetLoader(folder_path)
            % Konstruktor - ustawia główną ścieżkę do danych
            obj.FolderPath = folder_path;
            obj.Patients = [];
        end
        
        function load_data(obj)
            % Pobiera listę wszystkich elementów w folderze
            all_items = dir(obj.FolderPath);
            
            % Filtruje tylko katalogi i odrzuca systemowe '.' oraz '..'
            is_dir = [all_items.isdir];
            dir_names = {all_items.name};
            valid_dirs = is_dir & ~ismember(dir_names, {'.', '..'});
            patient_folders = all_items(valid_dirs);
            
            num_patients = length(patient_folders);
            
            if num_patients == 0
                warning('Nie znaleziono folderów pacjentów w ścieżce: %s', obj.FolderPath);
                return;
            end
            
            % Inicjalizacja pustej tablicy na obiekty
            obj.Patients = PatientData.empty(num_patients, 0);
            
            for i = 1:num_patients
                folder_name = patient_folders(i).name;
                patient_path = fullfile(obj.FolderPath, folder_name);
                
                % Ładowanie poszczególnych plików z folderu pacjenta
                % Zakładamy, że pliki nazywają się np. 'dep.mat', 'A12.mat' itd.

                dst = obj.load_file_or_empty(patient_path, 'ventricle.dst3d');
                dep = obj.load_file_or_empty(patient_path, 'user.dep');
                rep = obj.load_file_or_empty(patient_path, 'user.rep');
                A12 = obj.load_file_or_empty(patient_path, 'ventricles2standard_12.mat');
                A64 = obj.load_file_or_empty(patient_path, 'ventricles2BSM_(nijmegen_64).mat');
                y12 = obj.load_file_or_empty(patient_path, 'standard_12.refECG');
                y64 = obj.load_file_or_empty(patient_path, 'BSM_(nijmegen_64).refECG');
                
                % Utworzenie obiektu PatientData (jako ID używamy nazwy folderu)
                obj.Patients(i) = PatientData(folder_name, dep, rep, A12, A64, y12, y64);
            end
            
            fprintf('Pomyślnie załadowano dane %d pacjentów.\n', length(obj.Patients));
        end
        
        function patient = get_patient(obj, index_or_id)
            % Pobiera pacjenta po indeksie (liczba) lub po ID (string/char)
            if isnumeric(index_or_id)
                patient = obj.Patients(index_or_id);
            elseif ischar(index_or_id) || isstring(index_or_id)
                idx = strcmp({obj.Patients.PatientID}, index_or_id);
                patient = obj.Patients(idx);
                if isempty(patient)
                    error('Nie znaleziono pacjenta o ID: %s', index_or_id);
                end
            else
                error('Podaj poprawny indeks (liczba) lub ID (tekst).');
            end
        end
        
        function n = get_num_patients(obj)
            % Zwraca liczbę załadowanych pacjentów
            n = length(obj.Patients);
        end
    end
    
    methods (Access = private)
        function val = load_file_or_empty(~, folder_path, file_name)
            % Buduje pełną ścieżkę do pliku na podstawie podanej nazwy z rozszerzeniem
            file_path = fullfile(folder_path, file_name);
            
            % Sprawdza, czy plik istnieje
            if exist(file_path, 'file') == 2
                try
                    % Ładujemy plik (MATLAB radzi sobie z różnymi formatami przez 'load')
                    % Opcjonalnie można wymusić load(file_path, '-mat') jeśli pliki
                    % z niestandardowymi rozszerzeniami (.dep, .rep) są plikami MAT.
                    data = loadmat(file_path);
                    
                    if isstruct(data)
                        % Jeśli load zwrócił strukturę (np. plik .mat), 
                        % bierzemy pierwszą dostępną zmienną.
                        fields = fieldnames(data);
                        if ~isempty(fields)
                            val = data.(fields{1});
                        else
                            val = [];
                        end
                    else
                        % Jeśli load zwrócił bezpośrednio macierz (np. plik ASCII)
                        val = data;
                    end
                catch ME
                    % Zabezpieczenie na wypadek uszkodzonego pliku
                    warning('Plik istnieje, ale wystąpił błąd podczas ładowania: %s\n%s', ...
                            file_path, ME.message);
                    val = [];
                end
            else
                % Jeśli plik nie istnieje, zwracamy pustą tablicę
                val = [];
            end
        end
    end
end