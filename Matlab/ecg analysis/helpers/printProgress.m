function printProgress(iter, totalIter, suffix)
    percent = (iter / totalIter) * 100;

    barLength = 40;
    nFilled = floor(barLength * iter / totalIter);
    bar = [repmat('=', 1, nFilled), repmat(' ', 1, barLength - nFilled)];

    % Tekst do wypisania
    txt = sprintf('[%s] %6.2f%% (%d/%d) %s', bar, percent, iter, totalIter, suffix);

    % Liczba znaków tekstu – potrzebna do wymazania poprzedniego
    persistent prevLen;
    if isempty(prevLen)
        prevLen = 0;
    end

    % Wypisz odpowiednią liczbę backspace
    fprintf(repmat('\b', 1, prevLen));

    % Wypisz nowy pasek
    fprintf('%s', txt);

    % Zapamiętaj długość
    prevLen = length(txt);

    % Po zakończeniu przejdź do nowej linii i zresetuj
    if iter == totalIter
        fprintf('\n');
        prevLen = 0;
    end
end
