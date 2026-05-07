function meas = my_inversefunc(GEOM, initdep, initrep, mudep, murep)
% MY_INVERSEFUNC - Simplified inverse problem for depolarization and repolarization
% 
% INPUTS:
%   GEOM    - Wcześniej przygotowana struktura geometrii i sygnału (z loadera)
%   initdep - Ogniska/czasy początkowe dla depolaryzacji (np. z multifociscan)
%   initrep - Ogniska/czasy początkowe dla repolaryzacji
%   mudep   - Parametr regularyzacji (Laplasjan) dla depolaryzacji (np. 5e-6)
%   murep   - Parametr regularyzacji (Laplasjan) dla repolaryzacji (np. 1e-5)

disp('Inicjalizacja uproszczonego inversefunc...');

%% PARAMETRY GŁÓWNE
INV = struct();
INV.maxiter = 25;
INV.lambopt = 0.001;
INV.stopcrit = 2e-4;
INV.lambdamax = 500;
INV.MINRD = 0.15;
INV.mode = 4; % Mode 4 oznacza optymalizację Depolaryzacji i Repolaryzacji
INV.lpass = 1; % Zmiana z 5 na 1 (mniejsze wygładzanie, opcjonalnie)
INV.useWeighedRd = 0;

% Domyślne wartości z GEOM, bo ich potrzebujemy
INV.SPECS = GEOM.SPECS;
INV.AMA = GEOM.AMA;
INV.ATA = INV.AMA' * INV.AMA;

%% PRZYGOTOWANIE SYGNAŁU REFERENCYJNEGO
% Baseline correction and low-pass filtering
INV.BSM = baselinecor(lowpassma(GEOM.BSM(:, GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave), INV.lpass));
INV.PHIREF = INV.BSM;
INV.normphi = norm(INV.PHIREF, 'fro');

usetimes = size(INV.PHIREF, 2);
INV.T = ones(size(GEOM.AMA, 2), 1) * (0 : usetimes - 1);

%% REGULARYZACJA (Laplasjan na powierzchni siatki trójkątnej)
[INV.REGOP, INV.REGOPREP] = calcREGOP(GEOM, 1); % Używamy domyślnego surflapl
INV.ROTRO = INV.REGOP' * INV.REGOP;
INV.ROTROREP = INV.REGOPREP' * INV.REGOPREP;

%% INICJALIZACJA STRUKTUR OPTYMALIZACYJNYCH
OPT = struct();

% Struktura dla depolaryzacji
OPT.DEP.mu = mudep;
OPT.DEP.lambopt = INV.lambopt;
OPT.DEP.tims = initdep;
OPT.DEP.type = 'dep';

% Struktura dla repolaryzacji
OPT.REP.mu = murep;
OPT.REP.lambopt = INV.lambopt;
OPT.REP.tims = initrep;
OPT.REP.type = 'apd'; % 'apd' oznacza, że optymalizujemy odległość między dep a rep

% Puste struktury dla wyłączonych modułów (wymagane w argumentach funkcji pomocniczych)
OPT.NOT.pot = zeros(size(GEOM.VER,1),1); % Notch - wyłączony
OPT.AMP.pot = ones(size(GEOM.VER,1),1);  % Amplituda - wyłączona i ustawiona na 1

%% WSTĘPNA WERYFIKACJA STANU POCZĄTKOWEGO
TST = gettres_v(INV, OPT.DEP, OPT.REP, OPT.NOT, OPT.AMP);

disp('-------------------------------------------------------------------------');
disp(' iter    minD    maxD   stdD   minR    maxR   stdR     reg      rd     tres');
fprintf(' %3d %7.1f %7.1f %6.1f %6.1f %7.1f %6.1f %7.1f %7.4f %7.4f\n', ...
    0, min(OPT.DEP.tims), max(OPT.DEP.tims), std(OPT.DEP.tims), ...
    min(OPT.REP.tims), max(OPT.REP.tims), std(OPT.REP.tims), ...
    TST.reg, TST.rd, TST.tres);

%% GŁÓWNA PĘTLA ITERACYJNA MARQUARDTA
iter = 0;
depscore = 1;
repscore = 1;

while iter < INV.maxiter && (depscore || repscore) && TST.rd > INV.MINRD
    iter = iter + 1;
    
    % Krok 1: Optymalizacja Depolaryzacji
    [depscore, OPT.DEP, TST] = optimizeDepRep(INV, OPT.DEP, OPT.REP, OPT.NOT, OPT.AMP);
    opt = OPT.DEP;
    
    fprintf('D%3d %7.1f %7.1f %6.1f %6.1f %7.1f %6.1f %7.1f %7.4f %7.4f\n', ...
        iter, min(opt.tims), max(opt.tims), std(opt.tims), ...
        min(OPT.REP.tims), max(OPT.REP.tims), std(OPT.REP.tims), ...
        TST.reg, TST.rd, TST.tres);

    % Krok 2: Optymalizacja Repolaryzacji (jeśli jest potrzebna)
    [repscore, OPT.REP, TST] = optimizeDepRep(INV, OPT.REP, OPT.DEP, OPT.NOT, OPT.AMP);
    opt = OPT.REP;

    fprintf('R%3d %7.1f %7.1f %6.1f %6.1f %7.1f %6.1f %7.1f %7.4f %7.4f\n', ...
        iter, min(OPT.DEP.tims), max(OPT.DEP.tims), std(OPT.DEP.tims), ...
        min(opt.tims), max(opt.tims), std(opt.tims), ...
        TST.reg, TST.rd, TST.tres);
end

disp('-------------------------------------------------------------------------');

%% ZWRÓCENIE WYNIKÓW
% Ostatnia kontrola parametrów
TST = gettres_v(INV, OPT.DEP, OPT.REP, OPT.NOT, OPT.AMP);

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


function TST = my_gettres_v(INV, opt, keepopt, notchopt, amplopt)
    % Set weighted RD value parameter to 0.0010;
    wrd_param   = 0.0010;
    
    % Determine S for given .tims:
    % Check if type is depolarization or something else ('rep' or 'apd'):
    if strcmp(opt.type,'dep')
        tims_1 = opt.tims; tims_2 = keepopt.tims;
    else
        tims_1 = keepopt.tims; tims_2 = opt.tims;
    end
    
    % TST.S = getSfunc(INV.T, tims_1, tims_2, INV.SPECS, notchopt.pot, amplopt.pot, INV.mode, INV);   %
    TST.S = getSmode(INV.T, tims_1, tims_2, INV.SPECS, INV.mode);
    
    % Determine simulated BSM, residue, RD and weighted RD values:
    TST.PHIA    = lowpassma(INV.AMA * TST.S,INV.lpass);                                             % Construct filterd simulated BSM with A-matrix (INV.AMA) and simulated cardiac signals (TST.S)
    TST.RES     = INV.PHIREF-TST.PHIA(1:size(INV.PHIREF,1),1:size(INV.PHIREF,2));                   % Calculated residue --> difference INV.PHIREF and TST.PHIA
    TST.rd      = norm(TST.RES,'fro')/INV.normphi;                                                  % Calculate RD value
    TST.wrd     = sum(rms(TST.RES)./(wrd_param + rms(INV.PHIREF)));                                 % Calcualte weighted RD value
    
    % Determine regularization value:
    if strcmp(opt.type,'dep')
        TST.reg = norm(INV.REGOP*opt.tims);                                                         % Calculate regularization by multiplying REGOP with times
    else
        TST.reg = norm(INV.REGOPREP*opt.tims);                                                      % Calculate regularization by multiplying REGOPREP with times
    end
    
    % Determine value of regularization term plus RD value:
    if INV.useWeighedRd                                                                            % Check if weighted RD is used:
        if strcmp(opt.type,'dep')
            TST.tres = sqrt(TST.wrd^2+(TST.reg*opt.mu)^2);
        else
            TST.tres = sqrt(TST.wrd^2+(TST.reg*opt.mu)^2);
        end
    else
        if strcmp(opt.type,'dep')
            TST.tres = sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
        else
            TST.tres = sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
        end
    end
end