function TST = gettres_v_rep(INV, params_cell, fixed_dep, notchopt, amplopt)
    % GETTRES_V_1PARAM - Single parameter optimization with fixed depolarization
    
    wrd_param = 0.0010;
    
    % 1. Extract times (depolarization is static, repolarization is dynamic)
    dep_times = fixed_dep;
    rep_times = params_cell{1};
    
    % Safe defaults for optional arguments
    if nargin < 4 || isempty(notchopt)
        notchopt = struct('pot', zeros(length(dep_times), 1));
    end
    if nargin < 5 || isempty(amplopt)
        amplopt = struct('pot', ones(length(dep_times), 1));
    end
    
    % 2. Generate action potential matrix (S)
    TST.S = my_regional_Smode(INV, dep_times, rep_times);
    
    % 3. ECG simulation (PHIA) and signal residuals
    TST.PHIA = lowpassma(INV.AMA * TST.S, INV.lpass);
    TST.RES  = INV.PHIREF - TST.PHIA(1:size(INV.PHIREF,1), 1:size(INV.PHIREF,2));
    TST.rd   = norm(TST.RES, 'fro') / INV.normphi;
    TST.wrd  = sum(rms(TST.RES) ./ (wrd_param + rms(INV.PHIREF)));
    
    % 4. Regularization error calculation
    % Depolarization is fixed, thus we omit INV.REGOP multiplication to save compute time
    reg_rep = norm(INV.REGOPREP * rep_times);
    TST.reg = reg_rep; 
    
    % 5. Fetch regularization weights from INV structure
    mu_rep = 0;
    if isfield(INV, 'MuValues') && ~isempty(INV.MuValues)
        mu_rep = INV.MuValues(1); % Since we have 1 parameter, index is 1
    end
    
    % 6. Cost function calculation
    penalty_sq = (reg_rep * mu_rep)^2;
    
    if isfield(INV, 'useWeighedRd') && INV.useWeighedRd
        TST.tres = sqrt(TST.wrd^2 + penalty_sq);
    else
        TST.tres = sqrt(TST.rd^2 + penalty_sq);
    end
end