function meas = my_inversefunc(GEOM, initdep, initrep, mudep, murep)
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

%% WSTĘPNA WERYFIKACJA STANU POCZĄTKOWEGO
TST = gettres_v_from_TmpLut(INV, OPT.DEP, OPT.REP, OPT.NOT, OPT.AMP);

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
    if false
        % Visualize the current simulated ECG vs measured ECG
        figure(99); clf;
        set(gcf, 'Name', ['Optimization Progress - Iteration: ', num2str(iter)]);
        
        % Plot measured signal (Blue)
        sigplot(INV.PHIREF, 'Blue = Measured, Red = Simulated', GEOM.LAY, 1.5, 'b', 1);
        hold on;
        % Plot current simulated signal (Red)
        sigplot(TST.PHIA, '', GEOM.LAY, 1.5, 'r', 1);
        
        % Force MATLAB to update the figure window immediately
        drawnow; 
        % -------------------------------
    end
end

disp('-------------------------------------------------------------------------');

% Restore best model
disp(['Restoring best model from iteration ', num2str(best_iter), ' with RD = ', num2str(best_rd)]);
TST = gettres_v_from_TmpLut(INV, best_OPT_DEP, best_OPT_REP, OPT.NOT, OPT.AMP);

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