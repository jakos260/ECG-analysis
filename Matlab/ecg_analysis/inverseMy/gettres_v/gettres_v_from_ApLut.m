function TST = gettres_v_from_ApLut(INV, opt, keepopt, notchopt, amplopt)
    wrd_param = 0.0010;
    
    % 1. Rozróżnienie czasów (dep i rep)
    if strcmp(opt.type,'dep')
        tims_1 = opt.tims; 
        tims_2 = keepopt.tims;
    else
        tims_1 = keepopt.tims; 
        tims_2 = opt.tims;
    end
    
    % 2. Generowanie potencjałów przy użyciu naszego słownika LUT i regionów AHA
    TST.S = my_regional_Smode(INV, tims_1, tims_2);

    % 3. ORYGINALNA LOGIKA (BŁĘDY I REGULARYZACJA)
    TST.PHIA    = lowpassma(INV.AMA * TST.S, INV.lpass); 
    TST.RES     = INV.PHIREF - TST.PHIA(1:size(INV.PHIREF,1), 1:size(INV.PHIREF,2));
    TST.rd      = norm(TST.RES,'fro') / INV.normphi;
    TST.wrd     = sum(rms(TST.RES) ./ (wrd_param + rms(INV.PHIREF)));
    
    % Determine regularization value:
    if strcmp(opt.type,'dep')
        TST.reg = norm(INV.REGOP * opt.tims);
    else
        TST.reg = norm(INV.REGOPREP * opt.tims);
    end
    
    % Determine value of regularization term plus RD value:
    if INV.useWeighedRd
        TST.tres = sqrt(TST.wrd^2 + (TST.reg * opt.mu)^2);
    else
        TST.tres = sqrt(TST.rd^2 + (TST.reg * opt.mu)^2);
    end
end
