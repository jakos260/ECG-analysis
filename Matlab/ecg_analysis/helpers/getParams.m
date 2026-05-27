function params = getParams(valueSets, iter, totalIter, suffix)
% getParamsForIteration - wybiera zestaw parametrów dla danej iteracji
%
% valueSets - 1xN komórkowa tablica, każda komórka to wektor możliwych
%             wartości danego parametru (u Ciebie N = 5)
% iter      - numer iteracji (1..prod(liczby_wartości))
%
% params    - 1xN wektor wartości parametrów dla podanej iteracji

    % liczba parametrów (może być 5, ale funkcja działa ogólnie)
    nParams = numel(valueSets);

    % liczba możliwych wartości dla każdego parametru
    sizes = cellfun(@numel, valueSets);

    % całkowita liczba kombinacji
    totalComb = prod(sizes);

    if iter < 1 || iter > totalComb
        error('Numer iteracji musi być w zakresie od 1 do %d.', totalComb);
    end

    % konwersja iteracji (1-based) na "cyfry" w mieszanym systemie o podstawach = sizes
    idx0 = iter - 1;  % przechodzimy na indeks 0-based
    indices = zeros(1, nParams);

    % tutaj przyjmuję, że OSTATNI parametr zmienia się najszybciej
    for k = nParams:-1:1
        indices(k) = mod(idx0, sizes(k)) + 1;   % 1-based indeks dla parametru k
        idx0 = floor(idx0 / sizes(k));
    end

    % wybieramy wartości
    params = zeros(1, nParams);
    for k = 1:nParams
        params(k) = valueSets{k}(indices(k));
    end

    printProgress(iter, totalIter, suffix);
end


