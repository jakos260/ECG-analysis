classdef INVERSE_MODEL < handle
    % INVERSE_MODEL - A generalized inverse problem optimization class
    % utilizing my_optimize_nparam for N-parameter spatial optimization.
    
    properties
        Settings     % General solver settings
        ParamConfigs % Array of structs for parameters
        Optimizer    % Instance of my_optimize_nparam
        OptFunc      % Function handle for evaluation
        GEOM         % Geometry and spatial data
        INV          % Internal structure for optimization state
        BestIter     % Best iteration achieved
        BestRD       % Best relative distance (RD)
        BestParams   % Cell array of best parameters
    end
    
    methods
        function obj = INVERSE_MODEL(GEOM, param_configs, settings, opt_func)
            % Constructor
            disp('-------------------------------------------------------------------------');
            disp('Initializing generalized inverse_model...');
            obj.GEOM = GEOM;
            obj.ParamConfigs = param_configs;
            obj.Settings = settings;
            obj.OptFunc = opt_func;
            
            % Wypakowanie tablicy struktur (param_configs) do formatu cell
            num_params = length(param_configs);
            initial_params = cell(1, num_params);
            mu_values      = zeros(1, num_params);
            param_names    = cell(1, num_params);
            bounds         = cell(1, num_params);
            
            for i = 1:num_params
                initial_params{i} = param_configs(i).initial_value;
                mu_values(i)      = param_configs(i).mu;
                param_names{i}    = param_configs(i).name;
                bounds{i}         = param_configs(i).bounds;
            end
            
            % Konfiguracja struktury INV (bezpieczne, matlabowe przypisanie wartości domyślnych)
            obj.INV = struct();
            obj.INV.MuValues = mu_values;

            obj.INV.maxiter = 25;
            if isfield(settings, 'maxiter'); obj.INV.maxiter = settings.maxiter; end
            
            obj.INV.lambopt = 0.001;
            if isfield(settings, 'lambopt'); obj.INV.lambopt = settings.lambopt; end
            
            obj.INV.stopcrit = 2e-4;
            if isfield(settings, 'stopcrit'); obj.INV.stopcrit = settings.stopcrit; end
            
            obj.INV.MINRD = 0.15;
            if isfield(settings, 'MINRD'); obj.INV.MINRD = settings.MINRD; end
            
            obj.INV.lpass = 1;
            if isfield(settings, 'lpass'); obj.INV.lpass = settings.lpass; end
            
            obj.INV.mode = 4;
            if isfield(settings, 'mode'); obj.INV.mode = settings.mode; end
            
            obj.INV.useWeighedRd = 0;
            if isfield(settings, 'useWeighedRd'); obj.INV.useWeighedRd = settings.useWeighedRd; end
            
            obj.INV.lambdamax = 500;
            
            % Pozostałe parametry z geometrii
            obj.INV.SPECS = GEOM.SPECS;
            obj.INV.AMA = GEOM.AMA;
            obj.INV.ATA = obj.INV.AMA' * obj.INV.AMA;
            obj.INV.RegionIdx = GEOM.RegionIdx;
            
            if isfield(GEOM, 'LUT')
                obj.INV.LUT = GEOM.LUT;
            end
            
            % Przygotowanie sygnału referencyjnego i macierzy czasu T
            if isfield(GEOM, 'BSM') && isfield(GEOM.SPECS, 'onsetqrs') && isfield(GEOM.SPECS, 'endtwave')
                try
                    obj.INV.BSM = baselinecor(lowpassma(GEOM.BSM(:, GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave), obj.INV.lpass));
                catch
                     obj.INV.BSM = GEOM.BSM(:, GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave); % Fallback
                end
                obj.INV.PHIREF = obj.INV.BSM;
                obj.INV.normphi = norm(obj.INV.PHIREF, 'fro');
                
                usetimes = size(obj.INV.PHIREF, 2);
                obj.INV.T = ones(size(obj.INV.AMA, 2), 1) * (0 : usetimes - 1);
            else
                error('GEOM must contain BSM and SPECS.onsetqrs / endtwave');
            end
            
            % Inicjalizacja N-parametrowego optymalizatora
            obj.Optimizer = OPTIMIZE_NPARAM(...
                initial_params, mu_values, param_names, obj.INV.maxiter, bounds);
            
            % Mapowanie macierzy regularyzacji
            if isfield(GEOM, 'ROTRO_CELL')
                obj.INV.ROTRO_CELL = GEOM.ROTRO_CELL;
            else
                 try
                     % Zapisujemy REGOP i REGOPREP bezpośrednio do struktury INV,
                     % aby stara funkcja gettres_v_from_ApLut miała do nich dostęp.
                     [obj.INV.REGOP, obj.INV.REGOPREP] = calcREGOP(GEOM, 1);
                     obj.INV.ROTRO = obj.INV.REGOP' * obj.INV.REGOP;
                     obj.INV.ROTROREP = obj.INV.REGOPREP' * obj.INV.REGOPREP;
                 catch
                     disp('Warning: calcREGOP failed or not found, regularization might fail.');
                 end
            end
        end
        
        function meas = run_optimization(obj, notchopt, amplopt)
            if nargin < 2; notchopt = struct('pot', zeros(size(obj.GEOM.VER,1),1)); end
            if nargin < 3; amplopt = struct('pot', ones(size(obj.GEOM.VER,1),1)); end
            
            plot_en = isfield(obj.Settings, 'plot_en') && obj.Settings.plot_en;
            
            PLOT_STATE = [];
            if plot_en
                PLOT_STATE = obj.init_progress_plot();
            end
            
            % Startowa ewaluacja
            TST = obj.OptFunc(obj.INV, obj.Optimizer.Params, notchopt, amplopt);
            
            obj.BestRD = TST.rd;
            obj.BestParams = obj.Optimizer.Params;
            obj.BestIter = 0;
            
            if plot_en
                PLOT_STATE = obj.update_progress_plot(PLOT_STATE, 0, TST, obj.Optimizer.Params, obj.GEOM.ITRI);
                obj.plot_ecg_signals(0, obj.INV.PHIREF, TST.PHIA, obj.GEOM.LAY);
            end
            
            % Główna pętla
            for iter = 1:obj.INV.maxiter
                if TST.rd <= obj.INV.MINRD
                    disp(['Target RD (', num2str(obj.INV.MINRD), ') reached. Stopping.']);
                    break;
                end
                
                [score, TST] = obj.Optimizer.optimize(obj.INV, notchopt, amplopt);
                
                if TST.rd < obj.BestRD
                    obj.BestRD = TST.rd;
                    obj.BestParams = obj.Optimizer.Params;
                    obj.BestIter = iter;
                end
                
                param_means = cellfun(@mean, obj.Optimizer.Params);
                param_stds = cellfun(@std, obj.Optimizer.Params);
                msg = sprintf('ITER: %d | RD: %7.4f | Tres: %7.4f | DEP/REP Means: %4.0f %4.0f | DEP/REP Std: %4.0f %4.0f |', iter, TST.rd, TST.tres, param_means, param_stds);
                disp(msg);
                
                if plot_en
                    PLOT_STATE = obj.update_progress_plot(PLOT_STATE, iter, TST, obj.Optimizer.Params, obj.GEOM.ITRI);
                    obj.plot_ecg_signals(iter, obj.INV.PHIREF, TST.PHIA, obj.GEOM.LAY);
                end
                
                if score == 0
                    disp('Optimization stalled, no parameters improved.');
                    break; 
                end
            end
            
            disp('-------------------------------------------------------------------------');
            disp(['Restoring best model from iteration ', num2str(obj.BestIter), ' with RD = ', num2str(obj.BestRD)]);
            
            obj.Optimizer.Params = obj.BestParams;
            TST = obj.OptFunc(obj.INV, obj.BestParams, notchopt, amplopt);
            
            meas = struct();
            meas.final_params = obj.BestParams;
            meas.rdfinal = TST.rd;
            meas.regfinal = TST.reg;
            meas.tresfinal = TST.tres;
            
            COR = corrcoef(TST.PHIA, obj.INV.PHIREF);
            meas.corfinal = COR(2,1);
            meas.iterfinal = obj.BestIter;
            meas.simulated_ecg = TST.PHIA;
            
            disp(['Optimization completed. Best achieved RD = ', num2str(meas.rdfinal)]);
        end
        
        function result = predict(obj, new_params, notchopt, amplopt)
             if nargin < 3; notchopt = struct('pot', zeros(size(obj.GEOM.VER,1),1)); end
             if nargin < 4; amplopt = struct('pot', ones(size(obj.GEOM.VER,1),1)); end
             
             disp('Running prediction (forward model) on provided parameters...');
             TST = obj.OptFunc(obj.INV, new_params, notchopt, amplopt);
             
             result = struct();
             result.simulated_ecg = TST.PHIA;
             result.RD = TST.rd;
             if isfield(TST, 'S')
                 result.S = TST.S; 
             else
                 result.S = [];
             end
        end
    end
    
    %% PRIVATE PLOTTING METHODS
    methods (Access = private)
        function PLOT_STATE = init_progress_plot(obj)
            maxiter = obj.INV.maxiter;
            numP = obj.Optimizer.NumParams;
            pNames = obj.Optimizer.ParamNames;
            
            fig_debug = figure(22); clf;
            set(fig_debug, 'Name', 'Optimization progress N-Param', 'Color', 'w');
            
            PLOT_STATE = struct();
            PLOT_STATE.hist_iter = 0:maxiter;
            PLOT_STATE.hist_rd   = NaN(1, maxiter + 1);
            
            PLOT_STATE.hist_mean = cell(1, numP);
            PLOT_STATE.hist_std = cell(1, numP);
            for p=1:numP
                PLOT_STATE.hist_mean{p} = NaN(1, maxiter + 1);
                PLOT_STATE.hist_std{p} = NaN(1, maxiter + 1);
            end
            
            subplot(2,1,1);
            PLOT_STATE.h_rd = plot(PLOT_STATE.hist_iter, PLOT_STATE.hist_rd, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
            title('Relative Distance (RD)'); xlim([0, maxiter]); grid on;
            
            subplot(2,1,2);
            hold on;
            colors = lines(numP);
            PLOT_STATE.h_means = cell(1, numP);
            PLOT_STATE.h_stds = cell(1, numP);
            for p=1:numP
                PLOT_STATE.h_means{p} = plot(PLOT_STATE.hist_iter, PLOT_STATE.hist_mean{p}, ...
                    '-o', 'Color', colors(p,:), 'LineWidth', 1.5, 'MarkerFaceColor', colors(p,:), ...
                    'DisplayName', ['mean ', pNames{p}]);
                PLOT_STATE.h_stds{p} = plot(PLOT_STATE.hist_iter, PLOT_STATE.hist_mean{p}, ...
                    '-*', 'Color', colors(p,:), 'LineWidth', 1.5, 'MarkerFaceColor', colors(p,:), ...
                    'DisplayName', ['std ', pNames{p}]);
            end
            title('Parameter statistic values across network');
            xlabel('Iteration');
            legend('show');
            xlim([0, maxiter]); grid on;
            hold off;
        end
        
        function PLOT_STATE = update_progress_plot(obj, PLOT_STATE, iter, TST, current_params, ~)
            idx = iter + 1; 
            PLOT_STATE.hist_rd(idx) = TST.rd;
            set(PLOT_STATE.h_rd, 'YData', PLOT_STATE.hist_rd);
            
            for p=1:obj.Optimizer.NumParams
                PLOT_STATE.hist_mean{p}(idx) = mean(current_params{p});
                PLOT_STATE.hist_std{p}(idx) = std(current_params{p});
                set(PLOT_STATE.h_means{p}, 'YData', PLOT_STATE.hist_mean{p});
                set(PLOT_STATE.h_stds{p}, 'YData', PLOT_STATE.hist_std{p});
            end
            drawnow;
        end
        
        function plot_ecg_signals(~, iter, PHIREF, PHIA, LAY)
            figure(21); clf;
            set(gcf, 'Name', ['Optimization results - Iteration: ', num2str(iter)]);
            try
                sigplot(PHIREF, 'Blue = Measured, Red = Simulated', LAY, 1.5, 'b', 1);
                hold on;
                sigplot(PHIA, '', LAY, 1.5, 'r', 1);
            catch
                plot(PHIREF(:, 1), 'b', 'LineWidth', 1.5); hold on;
                plot(PHIA(:, 1), 'r', 'LineWidth', 1.5);
                title('Blue = Measured, Red = Simulated (Lead 1 Fallback)');
            end
            drawnow;
        end
    end
end