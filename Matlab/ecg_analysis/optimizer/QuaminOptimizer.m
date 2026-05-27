classdef QuaminOptimizer < handle
    properties
        ModelFunc       % Handle to the simulation model: @(p, dep, rep, L, A)
        InitialGuess    % Initial parameter guess vector
        LowerBounds     % Lower bounds for parameters
        UpperBounds     % Upper bounds for parameters
        MaxIter = 1000  % Maximum LM iterations
        Tolerance = 1e-5% Convergence tolerance
    end
    
    methods
        function obj = QuaminOptimizer(model_func, initial_guess, lb, ub)
            % Constructor
            obj.ModelFunc = model_func;
            obj.InitialGuess = initial_guess(:).'; 
            
            % Optional lower bounds
            if nargin > 2 && ~isempty(lb)
                obj.LowerBounds = lb(:).';
            else
                obj.LowerBounds = [];
            end
            
            % Optional upper bounds
            if nargin > 3 && ~isempty(ub)
                obj.UpperBounds = ub(:).';
            else
                obj.UpperBounds = [];
            end
        end
        
        function [theta_opt, final_error, iter] = run(obj, dataset_loader, modality_type)
            % Runs global optimization across the entire dataset (all patients)
            
            patients = dataset_loader.Patients;
            num_patients = length(patients);
            
            if num_patients == 0
                error('Dataset loader contains no patients.');
            end
            
            fprintf('Starting global optimization over %d patients for %s modality...\n', ...
                    num_patients, modality_type);
            
            % Create an objective function closure that captures 'patients' and 'modality_type'
            obj_func = @(p, need_G) obj.calculate_global_residuals(p, need_G, patients, modality_type);
            
            % Execute the standalone Levenberg-Marquardt function with constraints
            [theta_opt, final_error, iter] = levenberg_marquardt(...
                obj_func, obj.InitialGuess, obj.LowerBounds, obj.UpperBounds, ...
                obj.MaxIter, obj.Tolerance);
            
            fprintf('Optimization finished successfully.\n');
            fprintf('Total Iterations: %d\n', iter);
            fprintf('Global Error Norm: %.4f\n', final_error);
        end
    end
    
    methods (Access = private)
        function [res_total, G_total] = calculate_global_residuals(obj, p, need_G, patients, modality_type)
            % Calculates concatenated residuals and Jacobian for all patients (one full epoch)
            res_total = [];
            G_total = [];
            
            for i = 1:length(patients)
                patient = patients(i);
                [A, y_ref] = patient.get_modality_data(modality_type);
                
                % Skip patient if data is missing
                if isempty(A) || isempty(y_ref)
                    continue;
                end
                
                % 1. Calculate the estimated signal
                y_est = obj.ModelFunc(p, patient.dep, patient.rep, patient.L, A);
                
                % 2. Calculate local residuals (flattened)
                res_i = y_ref(:) - y_est(:);
                res_total = [res_total; res_i]; %#ok<AGROW> 
                
                % 3. Calculate local Jacobian if requested
                if need_G
                    G_i = obj.calculate_local_jacobian(p, patient, A);
                    G_total = [G_total; G_i]; %#ok<AGROW> 
                end
            end
            
            if isempty(res_total)
                error('No valid data found across all patients for modality: %s', modality_type);
            end
        end
        
        function G_i = calculate_local_jacobian(obj, p, patient, A)
            % Numerically computes the Jacobian matrix for a single patient
            eps_val = 1e-3;
            npar = numel(p);
            N = size(A, 1) * patient.L; 
            G_i = zeros(N, npar);
            
            for k = 1:npar
                dp = zeros(size(p));
                dp(k) = eps_val;
                
                y_plus  = obj.ModelFunc(p + dp, patient.dep, patient.rep, patient.L, A);
                y_minus = obj.ModelFunc(p - dp, patient.dep, patient.rep, patient.L, A);
                
                G_i(:, k) = (y_plus(:) - y_minus(:)) / (2 * eps_val);
            end
        end
    end
end