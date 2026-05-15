function meas = my_inversefunc(GEOM, initdep, initrep, mudep, murep)
% MY_INVERSEFUNC - Simplified inverse problem for depolarization and repolarization
disp('Inicjalizacja uproszczonego inversefunc...');

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
    
    % Krok 1: Optymalizacja Depolaryzacji
    [depscore, OPT.DEP, TST] = my_optimizeDepRep(INV, OPT.DEP, OPT.REP, OPT.NOT, OPT.AMP);
    opt = OPT.DEP;
    
    fprintf('D%3d %7.1f %7.1f %6.1f %6.1f %7.1f %6.1f %7.1f %7.4f %7.4f\n', ...
        iter, min(opt.tims), max(opt.tims), std(opt.tims), ...
        min(OPT.REP.tims), max(OPT.REP.tims), std(OPT.REP.tims), ...
        TST.reg, TST.rd, TST.tres);

    % Krok 2: Optymalizacja Repolaryzacji 
    [repscore, OPT.REP, TST] = my_optimizeDepRep(INV, OPT.REP, OPT.DEP, OPT.NOT, OPT.AMP);
    opt = OPT.REP;

    fprintf('R%3d %7.1f %7.1f %6.1f %6.1f %7.1f %6.1f %7.1f %7.4f %7.4f\n', ...
        iter, min(OPT.DEP.tims), max(OPT.DEP.tims), std(OPT.DEP.tims), ...
        min(opt.tims), max(opt.tims), std(opt.tims), ...
        TST.reg, TST.rd, TST.tres);
end

disp('-------------------------------------------------------------------------');

% Ostatnia kontrola parametrów
TST = gettres_v_from_TmpLut(INV, OPT.DEP, OPT.REP, OPT.NOT, OPT.AMP);

meas = struct();
meas.depfinal = OPT.DEP.tims;
meas.repfinal = OPT.REP.tims;
meas.rdfinal = TST.rd;
meas.regfinal = TST.reg;
meas.tresfinal = TST.tres;
COR = corrcoef(TST.PHIA, INV.PHIREF);
meas.corfinal = COR(2,1);
meas.iterfinal = iter;

disp(['Zakończono optymalizację po ', num2str(iter), ' iteracjach. RD = ', num2str(meas.rdfinal)]);
end