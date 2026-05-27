function padded_signal = pad_signal_custom(signal, idx_start, idx_end, pad_left_len, target_len)
    % Sprawdzenie czy sygnał to wektor kolumnowy (dla ułatwienia konkatenacji robimy z niego wiersz)
    is_column = iscolumn(signal);
    if is_column
        signal = signal';
    end

    % Pobranie wartości odniesienia do paddingu
    val_left = signal(idx_start);
    % end - idx_end + 1 oznacza n-tą próbkę od końca (np. dla 1 to jest dokładnie end)
    val_right = signal(end - idx_end + 1); 

    % Wygenerowanie paddingu z lewej strony
    pad_left = val_left * ones(1, pad_left_len);

    % Obliczenie ile próbek brakuje do docelowej długości
    % (zakładamy, że w środku zachowujemy cały oryginalny sygnał)
    current_total_len = length(pad_left) + length(signal);
    pad_right_len = target_len - current_total_len;

    % Zabezpieczenie przed podaniem zbyt małej docelowej długości
    if pad_right_len < 0
        error('Docelowa długość (%d) jest za mała! Sygnał + lewy padding zajmują już %d próbek.', target_len, current_total_len);
    end

    % Wygenerowanie paddingu z prawej strony
    pad_right = val_right * ones(1, pad_right_len);

    % Sklejenie wszystkiego w jeden wektor
    padded_signal = [pad_left, signal, pad_right];

    % Przywrócenie oryginalnego układu (kolumna/wiersz), jeśli to konieczne
    if is_column
        padded_signal = padded_signal';
    end
end