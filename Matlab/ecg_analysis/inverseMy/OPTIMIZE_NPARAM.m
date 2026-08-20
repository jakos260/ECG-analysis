classdef OPTIMIZE_NPARAM < handle
    % OPTIMIZE_NPARAM - A generalized N-parameter optimization class
    % based on Levenberg-Marquardt with spatial Tikhonov regularization.
    
    properties
        NumParams       % Number of parameters to optimize
        Params          % Cell array (1xN) of current parameter vectors
        MuValues        % Array (1xN) of regularization weights
        ParamNames      % Cell array (1xN) of parameter names
        MaxIter         % Maximum number of optimization iterations
        Bounds          % Cell array (1xN) of [min, max] bounds
        Lambdas         % Array (1xN) of current LM damping factors
    end
    
    methods
        function obj = OPTIMIZE_NPARAM(initial_params, mu_values, param_names, maxiter, bounds)
            % Constructor
            % Inputs:
            %   initial_params - Cell array of vectors containing initial times
            %   mu_values      - Array of mu weights for each parameter
            %   param_names    - (Optional) Cell array of strings for parameter names
            %   maxiter        - (Optional) Maximum global iterations (default = 25)
            %   bounds         - (Optional) Cell array of [min, max] for physiological clamping
            
            if ~iscell(initial_params)
                error('initial_params must be a cell array of parameter vectors.');
            end
            
            obj.NumParams = length(initial_params);
            obj.Params = initial_params;
            obj.MuValues = mu_values;
            
            % Handle optional param_names
            if nargin < 3 || isempty(param_names)
                obj.ParamNames = cell(1, obj.NumParams);
                for i = 1:obj.NumParams
                    obj.ParamNames{i} = sprintf('param_%d', i);
                end
            else
                if length(param_names) ~= obj.NumParams
                    error('Length of param_names must match initial_params.');
                end
                obj.ParamNames = param_names;
            end
            
            % Handle optional maxiter
            if nargin < 4 || isempty(maxiter)
                obj.MaxIter = 25;
            else
                obj.MaxIter = maxiter;
            end
            
            % Handle optional bounds (for Physiological Clamping)
            if nargin < 5 || isempty(bounds)
                obj.Bounds = cell(1, obj.NumParams);
                for i = 1:obj.NumParams
                    obj.Bounds{i} = [-Inf, Inf]; % Default: no bounds
                end
            else
                obj.Bounds = bounds;
            end
            
            % Initialize lambdas
            obj.Lambdas = ones(1, obj.NumParams) * 0.1; 
        end
        
        function [score, TST] = optimize(obj, INV, notchopt, amplopt)
            % OPTIMIZE - Sequentially optimizes all N parameters.
            % 
            % Note: This method calls generic wrapper functions:
            %  - gettres_v_Nparam(INV, params_cell, notchopt, amplopt)
            %  - my_regional_Smode_Nparam(INV, params_cell)
            % You will need to implement these wrappers in your project to map
            % the generic cell array (obj.Params) back to your specific structural 
            % layout (e.g., separating it back into opt, keepopt, etc.).
            
            score = 0;
            dif_tims = 1;
            
            % Initial evaluation
            TST = gettres_v_nparams(INV, obj.Params, notchopt, amplopt);
            
            for iter = 1:obj.MaxIter
                improvement_in_iter = false;
                
                for p = 1:obj.NumParams
                    starttres = TST.tres;
                    lamb = obj.Lambdas(p);
                    
                    % 1. Determine times for derivative calculation
                    params_plus = obj.Params;
                    params_min  = obj.Params;
                    
                    params_plus{p} = params_plus{p} + dif_tims;
                    params_min{p}  = params_min{p} - dif_tims;
                    
                    % 2. Calculate temporal derivative (Jacobian approximation)
                    Splus = my_regional_Smode(INV, params_plus{:});
                    Smin  = my_regional_Smode(INV, params_min{:});
                    
                    Sprime = (Splus - Smin) / (dif_tims * 2);
                    SST    = Sprime * Sprime';
                    
                    % 3. Compute regularization matrix
                    % Attempt to dynamically find the correct regularization matrix in INV
                    ROTRO_IN = [];
                    rotro_field_name = ['ROTRO_', upper(obj.ParamNames{p})];
                    if isfield(INV, rotro_field_name)
                        ROTRO_IN = INV.(rotro_field_name);
                    elseif isfield(INV, 'ROTRO_CELL') && length(INV.ROTRO_CELL) >= p
                        ROTRO_IN = INV.ROTRO_CELL{p};
                    elseif isfield(INV, 'ROTRO')
                        ROTRO_IN = INV.ROTRO; % Fallback
                    end
                    
                    if isempty(ROTRO_IN)
                        error('Could not find appropriate ROTRO regularization matrix in INV.');
                    end
                    
                    muREG = bsxfun(@times, obj.MuValues(p)^2, ROTRO_IN);
                    
                    % Hessian and Gradient
                    H = INV.ATA .* SST; 
                    GTGM = H + muREG; 
                    gtres = sum(INV.AMA' .* (Sprime * TST.RES'), 2); 
                    righths = gtres - muREG * obj.Params{p}; 
                    
                    stuck = 0;
                    maxstuck = 10;
                    param_updated = false;
                    
                    % 4. Iterate to find optimal LM damping parameter
                    while 1
                        lamb = lamb * 2;
                        
                        if lamb > INV.lambdamax
                            disp(['Optimization break: lambda reached maximum limit, lambda=', num2str(lamb)]); 
                            break; 
                        end
                        
                        % Damping factor for Levenberg-Marquardt
                        GTGM_L = GTGM + lamb^2 * eye(size(INV.AMA, 2));
                        
                        % Solve for step direction
                        deltau = GTGM_L \ righths;
                        
                        newtime = obj.Params{p} + deltau;
                        
                        % =========================================================================
                        % PHYSIOLOGICAL CLAMPING (BOUNDS)
                        % =========================================================================
                        min_b = obj.Bounds{p}(1);
                        max_b = obj.Bounds{p}(2);
                        newtime = max(min_b, min(newtime, max_b));
                        % =========================================================================
                        
                        test_params = obj.Params;
                        test_params{p} = newtime;
                        
                        test_TST = gettres_v_nparams(INV, test_params, notchopt, amplopt);
                        
                        if test_TST.tres < starttres
                            if (starttres - test_TST.tres) / starttres >= INV.stopcrit
                                obj.Params{p} = newtime;
                                TST = test_TST;
                                param_updated = true;
                                improvement_in_iter = true;
                                break;
                            else
                                stuck = stuck + 1;
                                if stuck >= maxstuck
                                    disp(['Optimization break: stuck parameter reached ', num2str(stuck)]); 
                                    break; 
                                end
                            end
                        else
                            stuck = 0;
                        end
                    end
                    
                    if param_updated
                        score = 1;
                        obj.Lambdas(p) = lamb / 4; % Update lambda for next cycle
                    end
                end
                
                % Early stopping if no parameter was improved in this full cycle
                if ~improvement_in_iter
                    break;
                end
            end
        end
    end
end