function S = my_regional_Smode(INV, dep, rep)
    % INV - struktura zawierająca INV.T, INV.RegionIdx oraz INV.LUT
    % dep - wektor czasów depolaryzacji [1500 x 1]
    % rep - wektor czasów repolaryzacji [1500 x 1]

    num_nodes = size(INV.T, 1);
    num_samples = size(INV.T, 2);
    S = zeros(num_nodes, num_samples);
    t_eval = INV.T(1, :); % globalny wektor próbkowania EKG
    
    num_regions = max(INV.RegionIdx);
    all_apds = [INV.LUT.APD]; % Wszystkie pre-kalkulowane długości
    
    for r = 1:num_regions
        % Wyizoluj węzły z obecnego klastra (AHA)
        nodes_in_region = find(INV.RegionIdx == r);
        if isempty(nodes_in_region), continue; end
        
        % 1. Znajdź średnie oczekiwane APD dla całego klastra
        avg_dep = mean(dep(nodes_in_region));
        avg_rep = mean(rep(nodes_in_region));
        target_APD = avg_rep - avg_dep;
        if target_APD < 1, target_APD = 1; end
        
        % 2. Wybierz najbardziej fizjologiczny szablon z LUT dla tego klastra
        [~, best_idx] = min(abs(all_apds - target_APD));
        V_template = INV.LUT(best_idx).V;
        t_tt2 = INV.LUT(best_idx).t;
        t_dep_template = INV.LUT(best_idx).t_dep;
        APD_template = INV.LUT(best_idx).APD;
        
        % 3. Aplikacja idealnego kształtu do poszczególnych komórek
        for n = 1:length(nodes_in_region)
            node = nodes_in_region(n);
            dep_n = dep(node);
            rep_n = rep(node);
            APD_n = rep_n - dep_n;
            if APD_n < 1, APD_n = 1; end
            
            t_mapped = zeros(size(t_tt2));
            
            % Faza depolaryzacji (czyste przesunięcie, zachowuje upstroke)
            mask_dep = t_tt2 <= t_dep_template;
            t_mapped(mask_dep) = t_tt2(mask_dep) - t_dep_template + dep_n;
            
            % Faza repolaryzacji (mikro-rozciągnięcie na potrzeby gradientu)
            % Współczynnik (APD_n / APD_template) z reguły będzie wynosił np. 1.002, 
            % co pozwala Marquardtowi działać bez niszczenia biofizyki.
            mask_rep = t_tt2 > t_dep_template;
            t_mapped(mask_rep) = dep_n + (t_tt2(mask_rep) - t_dep_template) * (APD_n / APD_template);
            
            % Ekstrapolacja i zapis do wynikowej macierzy S [1500 x czas]
            S(node, :) = interp1(t_mapped, V_template, t_eval, 'linear', 0);
        end
    end
end