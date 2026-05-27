classdef DatasetLoader < handle
    properties
        FolderPath      % Path to the main folder containing patient subdirectories
        Patients        % Array of PatientData objects
    end
    
    methods
        function obj = DatasetLoader(folder_path)
            % Constructor - sets the main data path
            obj.FolderPath = folder_path;
            obj.Patients = [];
        end
        
        function load_data(obj)
            % Gets a list of all items in the folder
            all_items = dir(obj.FolderPath);
            
            % Filters only directories and excludes system '.' and '..'
            is_dir = [all_items.isdir];
            dir_names = {all_items.name};
            valid_dirs = is_dir & ~ismember(dir_names, {'.', '..'});
            patient_folders = all_items(valid_dirs);
            
            num_patients = length(patient_folders);
            
            if num_patients == 0
                warning('No patient folders found in path: %s', obj.FolderPath);
                return;
            end
            
            % Initialize an empty array for objects
            obj.Patients = PatientData.empty(num_patients, 0);
            
            for i = 1:num_patients
                folder_name = patient_folders(i).name;
                patient_path = fullfile(obj.FolderPath, folder_name);
                
                % Load individual files from the patient's folder
                % Using specific file names and extensions as mapped below
                params      = obj.load_params_from_json(patient_path, 'params.json');
                y12_full    = obj.load_file_or_empty(patient_path, 'standard_12.refECG');
                y64_full    = obj.load_file_or_empty(patient_path, 'BSM_(nijmegen_64).refECG');
                y12 = y12_full(:, params.offset:end);
                y64 = y64_full(:, params.offset:end);
                A12 = obj.load_file_or_empty(patient_path, 'ventricles2standard_12.mat');
                A64 = obj.load_file_or_empty(patient_path, 'ventricles2BSM_(nijmegen_64).mat');
                dep = obj.load_file_or_empty(patient_path, 'user.dep');
                rep = obj.load_file_or_empty(patient_path, 'user.rep');
                dst = obj.load_file_or_empty(patient_path, 'ventricle.dst3d');
                
                % Create PatientData object (using folder name as ID)
                obj.Patients(i) = PatientData(folder_name, dep, rep, A12, A64, y12, y64, params);
            end
            
            fprintf('Successfully loaded data for %d patients.\n', length(obj.Patients));
        end
        
        function patient = get_patient(obj, index_or_id)
            % Retrieves a patient by index (numeric) or by ID (string/char)
            if isnumeric(index_or_id)
                patient = obj.Patients(index_or_id);
            elseif ischar(index_or_id) || isstring(index_or_id)
                idx = strcmp({obj.Patients.PatientID}, index_or_id);
                patient = obj.Patients(idx);
                if isempty(patient)
                    error('Patient not found with ID: %s', index_or_id);
                end
            else
                error('Provide a valid index (numeric) or ID (text).');
            end
        end
        
        function n = get_num_patients(obj)
            % Returns the number of loaded patients
            n = length(obj.Patients);
        end
    end
    
    methods (Access = private)
        function val = load_file_or_empty(~, folder_path, file_name)
            % Builds the full file path based on the provided name with extension
            file_path = fullfile(folder_path, file_name);
            
            % Checks if the file exists
            if exist(file_path, 'file') == 2
                try
                    % Load the file (custom loadmat function handles different formats)
                    data = loadmat(file_path);
                    
                    if isstruct(data)
                        % If loaded data is a struct, take the first available variable
                        fields = fieldnames(data);
                        if ~isempty(fields)
                            val = data.(fields{1});
                        else
                            val = [];
                        end
                    else
                        % If loaded data is directly a matrix
                        val = data;
                    end
                catch ME
                    % Safeguard in case of a corrupted file
                    warning('File exists, but an error occurred during loading: %s\n%s', ...
                            file_path, ME.message);
                    val = [];
                end
            else
                % If the file does not exist, return an empty array
                val = [];
            end
        end

        function params = load_params_from_json(~, folder_path, file_name)
            file_path = fullfile(folder_path, file_name);
            if exist(file_path, 'file') == 2
                txt = fileread(file_path);
                params = jsondecode(txt);
                return
            end
            error(sprintf("Can not find %s", file_path));
        end
    end
end