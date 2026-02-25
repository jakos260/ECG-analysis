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
    end
    
    methods
        function obj = PatientData(id, dep, rep, A12, A64, y12, y64)
            % Konstruktor przypisujący dane
            obj.PatientID = id;
            obj.dep = dep(:); % Wymuszenie wektora kolumnowego
            obj.rep = rep(:);
            obj.A12 = A12;
            obj.A64 = A64;
            obj.y12 = y12;
            obj.y64 = y64;
            
            % Automatyczne wyznaczenie długości sygnału L 
            % (zakładamy układ: [liczba kanałów x liczba próbek])
            if ~isempty(y12)
                obj.L = size(y12, 2);
            elseif ~isempty(y64)
                obj.L = size(y64, 2);
            else
                obj.L = 0;
            end
            
            % Opcjonalnie: Tutaj można dodać asercje sprawdzające 
            % czy wymiary A12 i y12 (oraz A64 i y64) są ze sobą zgodne.
        end
        
        function [A, y_ref] = get_modality_data(obj, modality_type)
            % Metoda pomocnicza zwracająca zestaw macierzy pod konkretny typ
            % Użyteczne, by optymalizator nie musiał wiedzieć, co przetwarza.
            switch modality_type
                case '12-lead'
                    A = obj.A12;
                    y_ref = obj.y12;
                case '64-lead'
                    A = obj.A64;
                    y_ref = obj.y64;
                otherwise
                    error('Nieobsługiwany typ modalności. Wybierz "12-lead" lub "64-lead".');
            end
        end
    end
end