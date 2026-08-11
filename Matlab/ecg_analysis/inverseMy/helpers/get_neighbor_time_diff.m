function [mean_diff, max_diff, all_diffs] = get_neighbor_time_diff(times, itri)
    % GET_NEIGHBOR_TIME_DIFF Oblicza różnice czasów między połączonymi węzłami.
    % 
    % Wejście:
    %   times - wektor czasów aktywacji (np. 1253x1)
    %   itri  - macierz trójkątów (Nx3)
    %
    % Wyjście:
    %   mean_diff - średnia absolutna różnica czasów [ms]
    %   max_diff  - maksymalna absolutna różnica czasów [ms]
    %   all_diffs - wektor wszystkich unikalnych różnic do głębszej analizy

    % 1. Wyciągnięcie wszystkich par węzłów (krawędzi) z macierzy trójkątów.
    % Każdy trójkąt [n1, n2, n3] generuje trzy krawędzie.
    edges = [itri(:,1), itri(:,2); 
             itri(:,2), itri(:,3); 
             itri(:,3), itri(:,1)];
             
    % 2. Standaryzacja i usunięcie duplikatów.
    % Sortujemy każdy wiersz (aby krawędź A-B była traktowana tak samo jak B-A),
    % a następnie wybieramy unikalne pary.
    edges = sort(edges, 2);
    unique_edges = unique(edges, 'rows');
    
    % 3. Obliczenie bezwzględnej różnicy czasów dla każdej unikalnej krawędzi.
    % times(unique_edges(:,1)) to czasy węzłów "początkowych",
    % times(unique_edges(:,2)) to czasy węzłów "końcowych".
    all_diffs = abs(times(unique_edges(:,1)) - times(unique_edges(:,2)));
    
    % 4. Agregacja wyników
    mean_diff = mean(all_diffs);
    max_diff = max(all_diffs);
end