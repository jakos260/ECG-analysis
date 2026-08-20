function TST = gettres_v_nparams(INV, params_cell, notchopt, amplopt)
    % GETTRES_V_NPARAMS - Standalone forward evaluation function
    % Zastępuje gettres_v_from_ApLut dla optymalizatora opartego na cell arrays.
    
    wrd_param = 0.0010;
    
    % 1. Rozpakowanie czasów z tablicy cell (indeks 1: dep, indeks 2: rep)
    dep_times = params_cell{1};
    rep_times = params_cell{2};
    
    % Bezpieczne wartości domyślne dla argumentów opcjonalnych
    if nargin < 3 || isempty(notchopt)
        notchopt = struct('pot', zeros(length(dep_times), 1));
    end
    if nargin < 4 || isempty(amplopt)
        amplopt = struct('pot', ones(length(dep_times), 1));
    end
    
    % 2. Generowanie macierzy potencjałów czynnościowych (S)
    TST.S = my_regional_Smode(INV, dep_times, rep_times);
    
    % 3. Symulacja EKG (PHIA) i wyliczenie reszt sygnałowych
    TST.PHIA = lowpassma(INV.AMA * TST.S, INV.lpass);
    TST.RES  = INV.PHIREF - TST.PHIA(1:size(INV.PHIREF,1), 1:size(INV.PHIREF,2));
    TST.rd   = norm(TST.RES, 'fro') / INV.normphi;
    TST.wrd  = sum(rms(TST.RES) ./ (wrd_param + rms(INV.PHIREF)));
    
    % 4. Wyliczenie błędu regularyzacji (kary za brak gładkości przestrzennej)
    reg_dep = norm(INV.REGOP * dep_times);
    reg_rep = norm(INV.REGOPREP * rep_times);
    
    % Sumaryczny parametr 'reg' zapisywany wyłącznie do wyświetlania w logach
    TST.reg = reg_dep + reg_rep; 
    
    % 5. Pobranie wag regularyzacji (mu) ze struktury INV
    mu_dep = 0; mu_rep = 0;
    if isfield(INV, 'MuValues') && length(INV.MuValues) >= 2
        mu_dep = INV.MuValues(1);
        mu_rep = INV.MuValues(2);
    end
    
    % 6. Obliczenie funkcji celu (tres) - zawiera błąd sygnału + kary regularyzacyjne
    penalty_sq = (reg_dep * mu_dep)^2 + (reg_rep * mu_rep)^2;
    
    if isfield(INV, 'useWeighedRd') && INV.useWeighedRd
        TST.tres = sqrt(TST.wrd^2 + penalty_sq);
    else
        TST.tres = sqrt(TST.rd^2 + penalty_sq);
    end
end