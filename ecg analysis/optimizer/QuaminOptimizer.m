classdef QuaminOptimizer < handle
    properties
        ModelFunc       % Uchwyt do modelu: @(p, dep, rep, L, A)
        InitialGuess    % Wektor startowy parametrów
        ErrorFunc       % Zewnętrzna funkcja błędu, np. @(y, y_fit) rmse(y, y_fit)
        MaxIter = 1000  % Maksymalna liczba iteracji algorytmu LM
        Tolerance = 1e-5 % Tolerancja zmiany błędu względnego
        Noneg = 0       % Flaga: wymuszenie nieujemnych parametrów
    end
    
    methods
        function obj = QuaminOptimizer(model_func, initial_guess, error_func)
            obj.ModelFunc = model_func;
            obj.InitialGuess = initial_guess(:).'; % Zawsze wektor wierszowy
            
            % Konfiguracja funkcji błędu (domyślnie RMSE na wektorach kolumnowych)
            if nargin > 2 && ~isempty(error_func)
                obj.ErrorFunc = error_func;
            else
                obj.ErrorFunc = @(y_true, y_est) rmse(y_true(:), y_est(:));
            end
        end
        
        function [theta_opt, y_fit_matrix, final_error, iter] = run(obj, patient_data, modality_type)
            % 1. Pobranie danych pacjenta z obiektu
            [A, y_ref] = patient_data.get_modality_data(modality_type);
            dep = patient_data.dep;
            rep = patient_data.rep;
            L = patient_data.L;
            
            % Sprawdzenie czy pacjent ma dane dla tej modalności
            if isempty(y_ref) || isempty(A)
                error('Pacjent %s nie posiada kompletnych danych dla modalności: %s', ...
                      patient_data.PatientID, modality_type);
            end
            
            % 2. Spłaszczenie sygnału referencyjnego na wektor kolumnowy
            y = y_ref(:); 
            
            % 3. Utworzenie lokalnego wrappera modelu (tylko od parametru p)
            internal_model = @(p) obj.ModelFunc(p, dep, rep, L, A);
            
            parms = obj.InitialGuess;
            npar = length(parms);
            
            % --- POCZĄTEK ALGORYTMU LEVENBERGA-MARQUARDTA ---
            [yest, G] = obj.calculate_jacobian(internal_model, parms);
            res = y - yest;
            normresiter = norm(res);
            
            iter = 0;
            lambda = 0.001;
            
            while true
                iter = iter + 1;
                if iter > obj.MaxIter
                    break;
                end
                
                GTG = G' * G;
                gtres = G' * res;
                break1 = 0;
                
                while break1 == 0
                    lambda = 2 * lambda;
                    if lambda > 1e10
                        break1 = 1;
                        break;
                    end
                    
                    % Krok Levenberga-Marquardta
                    MAT = GTG + (lambda^2) * eye(npar);
                    delp = pinv(MAT) * gtres;
                    
                    if obj.Noneg == 1
                        testp = max(parms + delp', 0);
                    else
                        testp = parms + delp';
                    end
                    
                    % Test nowej wartości theta
                    yest_test = internal_model(testp);
                    yest_test = yest_test(:);
                    res_test = y - yest_test;
                    testnorm = norm(res_test);
                    
                    if testnorm < normresiter
                        break;
                    end
                end
                
                relres = abs(testnorm - normresiter) / normresiter;
                if relres < obj.Tolerance || break1 == 1
                    break;
                end
                
                parms = testp;
                normresiter = testnorm;
                lambda = lambda / 4;
                
                [~, G] = obj.calculate_jacobian(internal_model, parms);
            end
            % --- KONIEC ALGORYTMU ---
            
            % Zapisanie wyników
            theta_opt = parms;
            
            % Wygenerowanie ostatecznego dopasowanego sygnału 
            % Zwracamy w oryginalnym kształcie (np. 12 x L)
            y_fit_matrix = internal_model(theta_opt);
            
            % Obliczenie ostatecznego błędu (używając przypisanej funkcji błędu)
            final_error = obj.ErrorFunc(y_ref, y_fit_matrix);
        end
    end
    
    methods (Access = private)
        function [yest, G] = calculate_jacobian(~, internal_model, p)
            % Oblicza numerycznie Jakobian
            y_model = internal_model(p);
            yest = y_model(:); % Spłaszczamy na wektor (np. 12*L x 1)
            
            eps_val = 1e-3;
            npar = numel(p);
            N = numel(yest);
            G = zeros(N, npar);
            
            for k = 1:npar
                dp = zeros(size(p));
                dp(k) = eps_val;
                
                y_plus  = internal_model(p + dp);
                y_minus = internal_model(p - dp);
                
                % Wypełniamy k-tą kolumnę Jakobianu spłaszczoną różnicą
                G(:, k) = (y_plus(:) - y_minus(:)) / (2 * eps_val);
            end
        end
    end
end