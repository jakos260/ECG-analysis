function meas = my_inversefunc(GEOM, initdep, initrep, mudep, murep, plot_en)
    % MY_INVERSEFUNC - Simplified inverse problem for depolarization and repolarization
    
    disp('Initializing simplified inversefunc...');
    
    %% PARAMETRY GŁÓWNE
    INV = struct();
    INV.maxiter = 25;
    INV.lambopt = 0.001;
    INV.stopcrit = 2e-4;
    INV.lambdamax = 500;
    INV.MINRD = 0.15;
    INV.mode = 4; 
    INV.lpass = 1; 
    INV.useWeighedRd = 0;
    
    INV.SPECS = GEOM.SPECS;
    INV.AMA = GEOM.AMA;
    INV.ATA = INV.AMA' * INV.AMA;
    INV.RegionIdx = GEOM.RegionIdx;
    INV.LUT = GEOM.LUT;
    
    %% PRZYGOTOWANIE SYGNAŁU REFERENCYJNEGO
    INV.BSM = baselinecor(lowpassma(GEOM.BSM(:, GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave), INV.lpass));
    INV.PHIREF = INV.BSM;
    INV.normphi = norm(INV.PHIREF, 'fro');
    
    usetimes = size(INV.PHIREF, 2);
    INV.T = ones(size(GEOM.AMA, 2), 1) * (0 : usetimes - 1);
    
    %% REGULARYZACJA 
    [INV.REGOP, INV.REGOPREP] = calcREGOP(GEOM, 1); 
    INV.ROTRO = INV.REGOP' * INV.REGOP;
    INV.ROTROREP = INV.REGOPREP' * INV.REGOPREP;
    
    %% INICJALIZACJA WYKRESÓW
    PLOT_STATE = [];
    if plot_en
        PLOT_STATE = init_progress_plot(INV.maxiter);
    end
    
    %% INICJALIZACJA STRUKTUR OPTYMALIZACYJNYCH
    OPT = struct();
    OPT.DEP.mu = mudep; OPT.DEP.lambopt = INV.lambopt; OPT.DEP.tims = initdep; OPT.DEP.type = 'dep';
    OPT.REP.mu = murep; OPT.REP.lambopt = INV.lambopt; OPT.REP.tims = initrep; OPT.REP.type = 'apd';
    OPT.NOT.pot = zeros(size(GEOM.VER,1),1);
    OPT.AMP.pot = ones(size(GEOM.VER,1),1);
    
    % Checkpoint trackers
    best_rd = inf;
    best_OPT_DEP = OPT.DEP;
    best_OPT_REP = OPT.REP;
    best_iter = 0;
    
    %% WSTĘPNA WERYFIKACJA STANU POCZĄTKOWEGO (ITERACJA 0)
    TST = gettres_v_from_ApLut(INV, OPT.DEP, OPT.REP, OPT.NOT, OPT.AMP);
    
    if plot_en
        PLOT_STATE = update_progress_plot(PLOT_STATE, 0, TST, OPT.DEP, OPT.REP, GEOM.ITRI);
        plot_ecg_signals(0, INV.PHIREF, TST.PHIA, GEOM.LAY);
    end
    
    disp('-------------------------------------------------------------------------');
    disp(' iter    minD    maxD   stdD   minR    maxR   stdR     reg      rd     tres');
    fprintf(' %3d %7.1f %7.1f %6.1f %6.1f %7.1f %6.1f %7.1f %7.4f %7.4f\n', ...
        0, min(OPT.DEP.tims), max(OPT.DEP.tims), std(OPT.DEP.tims), ...
        min(OPT.REP.tims), max(OPT.REP.tims), std(OPT.REP.tims), ...
        TST.reg, TST.rd, TST.tres);
    
    iter = 0; depscore = 1; repscore = 1;
    
    %% GŁÓWNA PĘTLA ITERACYJNA MARQUARDTA
    while iter < INV.maxiter && (depscore || repscore) && TST.rd > INV.MINRD
        iter = iter + 1;
    
        % Step 1: Depolarization optimization
        [depscore, OPT.DEP, TST] = my_optimizeDepRep(INV, OPT.DEP, OPT.REP, OPT.NOT, OPT.AMP);
        opt = OPT.DEP;
    
        if TST.rd < best_rd
            best_rd = TST.rd;
            best_OPT_DEP = OPT.DEP;
            best_OPT_REP = OPT.REP;
            best_iter = iter;
        end
    
        fprintf('D%3d %7.1f %7.1f %6.1f %6.1f %7.1f %6.1f %7.1f %7.4f %7.4f\n', ...
            iter, min(opt.tims), max(opt.tims), std(opt.tims), ...
            min(OPT.REP.tims), max(OPT.REP.tims), std(OPT.REP.tims), ...
            TST.reg, TST.rd, TST.tres);
    
        % Step 2: Repolarization optimization
        [repscore, OPT.REP, TST] = my_optimizeDepRep(INV, OPT.REP, OPT.DEP, OPT.NOT, OPT.AMP);
        opt = OPT.REP;
    
        if TST.rd < best_rd
            best_rd = TST.rd;
            best_OPT_DEP = OPT.DEP;
            best_OPT_REP = OPT.REP;
            best_iter = iter;
        end
    
        fprintf('R%3d %7.1f %7.1f %6.1f %6.1f %7.1f %6.1f %7.1f %7.4f %7.4f\n', ...
            iter, min(OPT.DEP.tims), max(OPT.DEP.tims), std(OPT.DEP.tims), ...
            min(opt.tims), max(opt.tims), std(opt.tims), ...
            TST.reg, TST.rd, TST.tres);
    
        % --- REAL-TIME VISUALIZATION ---
        if plot_en
            PLOT_STATE = update_progress_plot(PLOT_STATE, iter, TST, OPT.DEP, OPT.REP, GEOM.ITRI);
            plot_ecg_signals(iter, INV.PHIREF, TST.PHIA, GEOM.LAY);
        end
    end
    
    disp('-------------------------------------------------------------------------');
    
    % Restore best model
    disp(['Restoring best model from iteration ', num2str(best_iter), ' with RD = ', num2str(best_rd)]);
    TST = gettres_v_from_ApLut(INV, best_OPT_DEP, best_OPT_REP, OPT.NOT, OPT.AMP);
    
    meas = struct();
    meas.depfinal = best_OPT_DEP.tims;
    meas.repfinal = best_OPT_REP.tims;
    meas.rdfinal = TST.rd;
    meas.regfinal = TST.reg;
    meas.tresfinal = TST.tres;
    COR = corrcoef(TST.PHIA, INV.PHIREF);
    meas.corfinal = COR(2,1);
    meas.iterfinal = best_iter;
    meas.simulated_ecg = TST.PHIA;
    
    disp(['Optimization completed. Best achieved RD = ', num2str(meas.rdfinal)]);

end

%% =======================================================================
%  PLOTTING FUNCTION
%  =======================================================================

function PLOT_STATE = init_progress_plot(maxiter)
    % INIT_PROGRESS_PLOT
    
    fig_debug = figure(22); clf;
    set(fig_debug, 'Name', 'Optimization progress', 'Color', 'w');
    
    PLOT_STATE = struct();
    PLOT_STATE.hist_iter     = 0:maxiter;
    PLOT_STATE.hist_rd       = NaN(1, maxiter + 1);
    PLOT_STATE.hist_mean_dep = NaN(1, maxiter + 1);
    PLOT_STATE.hist_max_dep  = NaN(1, maxiter + 1);
    PLOT_STATE.hist_mean_rep = NaN(1, maxiter + 1);
    PLOT_STATE.hist_max_rep  = NaN(1, maxiter + 1);
    
    subplot(5,1,1);
    PLOT_STATE.h_rd = plot(PLOT_STATE.hist_iter, PLOT_STATE.hist_rd, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    title('Relative Distance (RD)'); xlim([0, maxiter]); grid on;
    
    subplot(5,1,2);
    PLOT_STATE.h_mean_dep = plot(PLOT_STATE.hist_iter, PLOT_STATE.hist_mean_dep, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    title('Mean DEP Difference between Neighbors [ms]'); xlim([0, maxiter]); grid on;
    
    subplot(5,1,3);
    PLOT_STATE.h_max_dep = plot(PLOT_STATE.hist_iter, PLOT_STATE.hist_max_dep, 'b-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    title('Max DEP Difference between Neighbors [ms]'); xlim([0, maxiter]); grid on;
    
    subplot(5,1,4);
    PLOT_STATE.h_mean_rep = plot(PLOT_STATE.hist_iter, PLOT_STATE.hist_mean_rep, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    title('Mean REP Difference between Neighbors [ms]'); xlim([0, maxiter]); grid on;
    
    subplot(5,1,5);
    PLOT_STATE.h_max_rep = plot(PLOT_STATE.hist_iter, PLOT_STATE.hist_max_rep, 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    title('Max REP Difference between Neighbors [ms]'); xlim([0, maxiter]); grid on;
    xlabel('Iteration');
end

function PLOT_STATE = update_progress_plot(PLOT_STATE, iter, TST, OPT_DEP, OPT_REP, ITRI)
    % UPDATE_PROGRESS_PLOT
    
    [mean_d, max_d] = get_neighbor_time_diff(OPT_DEP.tims, ITRI);
    [mean_r, max_r] = get_neighbor_time_diff(OPT_REP.tims, ITRI);
    
    idx = iter + 1; % Zabezpieczenie dla iteracji 0
    PLOT_STATE.hist_rd(idx)       = TST.rd;
    PLOT_STATE.hist_mean_dep(idx) = mean_d;
    PLOT_STATE.hist_max_dep(idx)  = max_d;
    PLOT_STATE.hist_mean_rep(idx) = mean_r;
    PLOT_STATE.hist_max_rep(idx)  = max_r;
    
    set(PLOT_STATE.h_rd, 'YData', PLOT_STATE.hist_rd);
    set(PLOT_STATE.h_mean_dep, 'YData', PLOT_STATE.hist_mean_dep);
    set(PLOT_STATE.h_max_dep, 'YData', PLOT_STATE.hist_max_dep);
    set(PLOT_STATE.h_mean_rep, 'YData', PLOT_STATE.hist_mean_rep);
    set(PLOT_STATE.h_max_rep, 'YData', PLOT_STATE.hist_max_rep);
    drawnow;
end

function plot_ecg_signals(iter, PHIREF, PHIA, LAY)
    % PLOT_ECG_SIGNALS
    
    figure(21); clf;
    set(gcf, 'Name', ['Optimization results - Iteration: ', num2str(iter)]);
    
    sigplot(PHIREF, 'Blue = Measured, Red = Simulated', LAY, 1.5, 'b', 1);
    hold on;
    sigplot(PHIA, '', LAY, 1.5, 'r', 1);
    drawnow;

end