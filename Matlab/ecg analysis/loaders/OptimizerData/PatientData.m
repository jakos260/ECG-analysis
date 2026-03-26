classdef PatientData
    properties
        PatientID   % Unikalny identyfikator, np. 'Pat_001'
        dep         % Wektor czasów depolaryzacji węzłów [liczba_węzłów x 1]
        rep         % Wektor czasów repolaryzacji węzłów [liczba_węzłów x 1]
        A12         % Macierz transformacji dla 12 odprowadzeń [12 x liczba_węzłów]
        A64         % Macierz transformacji dla 64 odprowadzeń [64 x liczba_węzłów]
        y12         % Sygnał referencyjny 12-lead [12 x L]
        y64         % Sygnał referencyjny 64-lead [64 x L]
        L           % Długość sygnałów w próbkach (wyliczana automatycznie)
        params      % Offsety
    end
    
    methods
        function obj = PatientData(id, dep, rep, A12, A64, y12, y64, params)
            % Constructor assigning data
            obj.PatientID   = id;
            obj.params      = params;
            obj.dep = dep(:);
            obj.rep = rep(:);
            obj.A12 = A12;
            obj.A64 = A64;
            obj.y12 = y12;
            obj.y64 = y64;
            
            % Automatically determine signal length L 
            % (assuming layout: [num_channels x num_samples])
            if ~isempty(y12)
                obj.L = size(y12, 2);
            elseif ~isempty(y64)
                obj.L = size(y64, 2);
            else
                obj.L = 0;
            end
            
            % Optional: Add assertions here to check if dimensions of 
            % A12 and y12 (and A64/y64) match each other.
        end
        
        function [A, y_ref] = get_modality_data(obj, modality_type)
            % Helper method returning a set of matrices for a specific type.
            % Useful so the optimizer doesn't need to know what it's processing.
            switch modality_type
                case '12-lead'
                    A = obj.A12;
                    y_ref = obj.y12;
                case '64-lead'
                    A = obj.A64;
                    y_ref = obj.y64;
                otherwise
                    error('Unsupported modality type. Choose "12-lead" or "64-lead".');
            end
        end
    end
end