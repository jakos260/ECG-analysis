function TST = gettres_v_mapper(INV, params_cell, fixed_dep, notchopt, amplopt)
    % GETTRES_V_MAPPER - Wersja optymalizacji wykorzystująca gets_MapperVersion
    % (obliczanie potencjałów S zgodnie z formułą z C++ PCyncComputeECG::getS)
    
    wrd_param = 0.0010;
    
    % 1. Extract times (depolarization is static, repolarization is dynamic)[cite: 4]
    dep_times = fixed_dep;
    rep_times = params_cell{1};
    
    % Safe defaults for optional arguments[cite: 4]
    if nargin < 4 || isempty(notchopt)
        notchopt = struct('pot', zeros(length(dep_times), 1));
    end
    if nargin < 5 || isempty(amplopt)
        amplopt = struct('pot', ones(length(dep_times), 1));
    end
    
    % 2. Generate action potential matrix (S) using gets_MapperVersion
    maxt = size(INV.T, 2);
    TST.S = gets_mapper(dep_times, rep_times, maxt);
    
    % 3. ECG simulation (PHIA) and signal residuals[cite: 4]
    TST.PHIA = lowpassma(INV.AMA * TST.S, INV.lpass);
    TST.RES  = INV.PHIREF - TST.PHIA(1:size(INV.PHIREF,1), 1:size(INV.PHIREF,2));
    TST.rd   = norm(TST.RES, 'fro') / INV.normphi;
    TST.wrd  = sum(rms(TST.RES) ./ (wrd_param + rms(INV.PHIREF)));
    
    % 4. Regularization error calculation[cite: 4]
    % Depolarization is fixed, thus we omit INV.REGOP multiplication to save compute time
    reg_rep = norm(INV.REGOPREP * rep_times);
    TST.reg = reg_rep; 
    
    % 5. Fetch regularization weights from INV structure[cite: 4]
    mu_rep = 0;
    if isfield(INV, 'MuValues') && ~isempty(INV.MuValues)
        mu_rep = INV.MuValues(1); % Since we have 1 parameter, index is 1
    end
    
    % 6. Cost function calculation[cite: 4]
    penalty_sq = (reg_rep * mu_rep)^2;
    
    if isfield(INV, 'useWeighedRd') && INV.useWeighedRd
        TST.tres = sqrt(TST.wrd^2 + penalty_sq);
    else
        TST.tres = sqrt(TST.rd^2 + penalty_sq);
    end
end