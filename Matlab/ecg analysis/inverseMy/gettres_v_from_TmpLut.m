function TST = gettres_v_from_TmpLut(INV, opt, keepopt, notchopt, amplopt)
    wrd_param = 0.0010;
    
    % 1. Rozróżnienie czasów (dep i rep)
    if strcmp(opt.type,'dep')
        tims_dep = opt.tims; 
        tims_rep = keepopt.tims;
    else
        tims_dep = keepopt.tims; 
        tims_rep = opt.tims;
    end
    
    num_nodes = size(INV.T, 1);
    num_samples = size(INV.T, 2);
    TST.S = zeros(num_nodes, num_samples);
    t_eval = INV.T(1, :); % wektor czasu próbkowania EKG
    
    % Wymagane: INV.LUT (Twój pre-kalkulowany słownik)
    % Wymagane: INV.RegionIdx (Wektor 1500x1, przypisujący węzeł do ID klastra, np. 1-50)
    
    num_regions = max(INV.RegionIdx);
    
    % 2. APLIKACJA SŁOWNIKA DO KLASTRÓW
    for r = 1:num_regions
        % Znajdź wszystkie węzły należące do tego klastra
        nodes_in_region = find(INV.RegionIdx == r);
        
        % Oblicz średnie oczekiwane APD dla tego obszaru wg optymalizatora
        avg_dep = mean(tims_dep(nodes_in_region));
        avg_rep = mean(tims_rep(nodes_in_region));
        target_APD = avg_rep - avg_dep;
        
        % Znajdź w słowniku LUT przebieg o fizjologicznym APD najbliższym target_APD
        all_apds = [INV.LUT.APD];
        [~, best_lut_idx] = min(abs(all_apds - target_APD));
        
        % Pobieramy idealny, fizjologiczny kształt fazy 3
        V_template = INV.LUT(best_lut_idx).V;
        APD_template = INV.LUT(best_lut_idx).APD; % Czas trwania tego wzorca
        
        % Aplikujemy ten wzorzec do wszystkich węzłów w klastrze,
        % przesuwając go w czasie zgodnie z ICH WŁASNYM czasem depolaryzacji
        for n = 1:length(nodes_in_region)
            node = nodes_in_region(n);
            dep_n = tims_dep(node);
            
            % Mapowanie tylko w oparciu o przesunięcie (brak rozciągania fazy!)
            % Ponieważ V_template ma już wymaganą długość (wybraną ze słownika)
            t_mapped = t_eval - dep_n + (APD_template / 2); % Uproszczone wyrównanie
            
            % Tu wykorzystujemy Twoje dt=0.1s (wzorce ze słownika) do interpolacji na czas EKG
            t_tt2 = (0:length(V_template)-1) * 0.1; % czas oryginalnego generatora
            
            % Wpisanie do dużej macierzy 
            TST.S(node, :) = interp1(t_tt2 + dep_n, V_template, t_eval, 'linear', 0);
        end
    end

    % 3. ORYGINALNA LOGIKA (BŁĘDY I REGULARYZACJA)
    % (Dokładnie taka sama reszta jak w poprzedniej odpowiedzi...)
    TST.PHIA = lowpassma(INV.AMA * TST.S, INV.lpass); 
    % ... i tak dalej
end