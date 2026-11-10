function S = gets_mapper(dep, rep, maxt)
    % GETS_CPP - Zwraca macierz potencjałów S dla wektorów czasów depolaryzacji i repolaryzacji
    % na podstawie formuły C++ (PCyncComputeECG::getS)
    
    % --- ZHARDKODOWANE PARAMETRY KONFIGURACYJNE ---
    USEDEPSLOPE     = true; 
    DEPSLOPE        = -0.5;   % Zastąp docelową stałą nachylenia depolaryzacji
    
    mSourceVersion  = 1;     % Wybierz 1 lub 0 w zależności od wersji źródłowej repolaryzacji
    plateauSlope    = 0.2;   % Nachylenie plateau (odp. mBeat.getPlateauSlope())
    repSlope        = 0.3;   % Nachylenie repolaryzacji (odp. mBeat.getRepSlope())
    initialRepSlope = 0.0;   % Według komentarza w C++ ta wartość to zawsze 0
    
    % --- PRZYGOTOWANIE ZMIENNYCH ---
    nn = length(dep);
    nt = maxt;
    
    % Upewniamy się, że dep i rep to wektory kolumnowe[cite: 1]
    if size(dep, 1) ~= nn, dep = dep'; end
    if ~isempty(rep) && size(rep, 1) ~= nn, rep = rep'; end
    
    % Wektor czasu 0 do maxt-1 (odpowiednik pętli po "t" z C++)
    T = 0:(nt - 1);
    
    % Macierze różnicy czasu względem aktywacji/repolaryzacji węzłów[cite: 1]
    TDEP = T - dep * ones(1, nt); 
    
    % Inicjalizacja wynikowej macierzy potencjałów
    S = zeros(nn, nt);
    
    % Znajdujemy indeksy węzłów posiadających poprawne czasy (w C++: if dep[ir] >= 0)
    valid = (dep >= 0);
    
    % --- FAZA DEPOLARYZACJI ---
    if USEDEPSLOPE
        % Formuła z C++: rowS[t] = 1.0 /( 1.0 + exp(DEPSLOPE * tdep ) );
        S(valid, :) = 1.0 ./ (1.0 + exp(DEPSLOPE .* TDEP(valid, :)));
    else
        warning('USEDEPSLOPE ustawione na false - brak zewnętrznej funkcji computedSrowDep.');
    end
    
    % --- FAZA REPOLARYZACJI ---
    if ~isempty(rep)
        % Wyznaczamy czasy repolaryzacji macierzowo
        TREP = T - rep * ones(1, nt);
        
        if mSourceVersion == 1
            % Formuła C++ dla wersji 1
            vali = (1.0 ./ (1.0 + exp(plateauSlope .* TREP(valid, :)))) .* ...
                   (1.0 ./ (1.0 + exp(repSlope .* TREP(valid, :))));
            
            S(valid, :) = S(valid, :) .* vali;
            
        else % mSourceVersion == 0
            % Formuła C++ dla wersji 0 (wcześniej realizowana jako for wewnątrz for)
            y_inc = (initialRepSlope + ...
                    (1.0 ./ (1.0 + exp(plateauSlope .* TREP(valid, :))))) .* ...
                    (1.0 ./ (1.0 + exp(repSlope .* TREP(valid, :))));
            
            % Szybka kumulacja sumy wzdłuż wierszy (odp. sumY += ...)
            y = cumsum(y_inc, 2);
            sumY = y(:, end);
            
            % Skalowanie od 0 do 1
            S(valid, :) = S(valid, :) .* (1 - bsxfun(@rdivide, y, sumY));
        end
        
        % Normalizacja wyniku: podzielenie wierszy przez ich maksymalną amplitudę[cite: 1]
        maxAmpl = max(S(valid, :), [], 2);
        S(valid, :) = bsxfun(@rdivide, S(valid, :), maxAmpl);
    end
end