function TST = gettres_v_from_TT2mod(INV, opt, keepopt, notchopt, amplopt)
    % Set weighted RD value parameter to 0.0010;
    wrd_param = 0.0010;
    
    % Rozróżnienie wektorów czasowych w zależności od optymalizowanej fazy
    if strcmp(opt.type,'dep')
        tims_1 = opt.tims; 
        tims_2 = keepopt.tims;
    else
        tims_1 = keepopt.tims; 
        tims_2 = opt.tims;
    end
    
    % =========================================================================
    % 1. GENERACJA WZORCA AP (TEMPLATE) ZA POMOCĄ MODELU TEN TUSSCHER 2
    % =========================================================================
    % Parametry ustawione na sztywno według zaleceń:
    HT = 0.1; 
    STOPTIME = 500; 
    CT = 1; % 1 = epicardial (domyślny typ komórki)
    phases_mod = [1, 1, 1]; % Domyślne wagi prądów
    HR = 60; 
    ISO = 0;
    
    % Uruchamiamy generator TYLKO RAZ dla wydajności
    [t_tt2, V_tt2] = wrapper_TenTusscher2mod(HT, STOPTIME, CT, phases_mod, HR, ISO);
    
    % Ograniczenie potencjału do zakresu [0, 1] - macierz AMA (Transfer Matrix) 
    % w tej architekturze (EDL) zawsze oczekuje znormalizowanych dipoli (0 do 1)
    V_norm = (V_tt2 - min(V_tt2)) / (max(V_tt2) - min(V_tt2));
    
    % Wyciągamy punkty charakterystyczne wzorca (żeby wiedzieć jak go rozciągać)
    % Czas depolaryzacji wzorca (przekroczenie 50% amplitudy)
    idx_dep_temp = find(V_norm >= 0.5, 1, 'first');
    % Czas repolaryzacji wzorca (spadek poniżej 10% amplitudy - APD90)
    idx_rep_temp = find(V_norm(idx_dep_temp:end) <= 0.1, 1, 'first') + idx_dep_temp - 1;
    if isempty(idx_rep_temp), idx_rep_temp = length(V_norm); end
    
    t_dep_temp = t_tt2(idx_dep_temp);
    t_rep_temp = t_tt2(idx_rep_temp);
    APD_temp = t_rep_temp - t_dep_temp; % Czas trwania potencjału wzorca
    
    % =========================================================================
    % 2. MAPOWANIE WZORCA NA 1500 WĘZŁÓW SERCA (SHIFT & STRETCH)
    % =========================================================================
    num_nodes = size(INV.T, 1);
    num_samples = size(INV.T, 2);
    TST.S = zeros(num_nodes, num_samples); % Inicjalizacja pustej macierzy [1500x400]
    
    t_eval = INV.T(1, :); % Globalny wektor czasu EKG (np. od 0 do 400 ms)
    
    for i = 1:num_nodes
        dep_i = tims_1(i);
        rep_i = tims_2(i);
        APD_i = rep_i - dep_i;
        
        % Zabezpieczenie przed ujemnym/zerowym APD (koliduje przy dzieleniu)
        if APD_i < 1, APD_i = 1; end 
        
        t_mapped = zeros(size(t_tt2));
        
        % A) Faza przed depolaryzacją i upstroke (zwykłe przesunięcie w czasie)
        mask_dep = t_tt2 <= t_dep_temp;
        t_mapped(mask_dep) = t_tt2(mask_dep) - t_dep_temp + dep_i;
        
        % B) Faza plateau i repolaryzacja (przesunięcie oraz rozciągnięcie/skurczenie)
        mask_rep = t_tt2 > t_dep_temp;
        t_mapped(mask_rep) = dep_i + (t_tt2(mask_rep) - t_dep_temp) * (APD_i / APD_temp);
        
        % Interpolacja tak przekształconego wzorca na sztywne próbki czasu EKG
        % Ekstrapolacja '0' na zewnątrz (potencjał spoczynkowy)
        TST.S(i, :) = interp1(t_mapped, V_norm, t_eval, 'linear', 0);
    end

    % =========================================================================
    % 3. ORYGINALNA LOGIKA (BŁĘDY I REGULARYZACJA)
    % =========================================================================
    % Determine simulated BSM, residue, RD and weighted RD values:
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