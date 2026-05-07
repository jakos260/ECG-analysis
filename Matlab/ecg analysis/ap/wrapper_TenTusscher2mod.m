function [t, V] = wrapper_TenTusscher2mod(HT, STOPTIME, CT, phases_mod, HR, ISO)
% WRAPPER_TENTUSSCHER2MOD Wrapper dla modelu Ten Tusscher 2
% Generuje kilka uderzeń dla zadanego HR, aby zasymulować zjawisko 
% adaptacji komórki (restitution), a zwraca tylko OSTATNI potencjał.

    % --- Obsługa opcjonalnych argumentów ---
    if nargin < 5 || isempty(HR)
        HR = 60;
    end

    if nargin < 6 || isempty(ISO)
        ISO = 0;
    end

    % --- Parametry stymulacji ---
    % Liczba uderzeń (5 to świetny kompromis między dokładnością 
    % fizjologiczną a szybkością responsywności UI)
    N_beats = 5; 
    
    % Długość cyklu serca w milisekundach (BCL)
    BCL = 60000 / HR; 
    
    % Amplituda i czas trwania impulsów
    Stim_I = ones(1, N_beats) * -52;   
    Stim_T = ones(1, N_beats) * 1;     
    
    % Interwały: pierwszy impuls po 1 ms, kolejne co zadany BCL
    Stim_Int = [1, ones(1, N_beats - 1) * BCL]; 
    
    % Całkowity czas symulacji:
    % Musi pomieścić N-1 interwałów BCL oraz czas (STOPTIME) dla ostatniego uderzenia.
    % Dodajemy +2 ms marginesu bezpieczeństwa na przesunięcia początkowe.
    TOTAL_STOPTIME = (N_beats - 1) * BCL + STOPTIME + 2;

    % --- Wywołanie głównej funkcji modelu ---
    [t_full, V_full] = TenTusscher2mod(HT, TOTAL_STOPTIME, Stim_I, Stim_T, ISO, Stim_Int, CT, phases_mod);

    % --- Wycięcie ostatniego uderzenia ---
    % Bazowy czas, w którym rozpoczęło się ostatnie uderzenie
    last_beat_start = (N_beats - 1) * BCL;
    
    % Znajdujemy indeksy odpowiadające ostatniemu przebiegowi
    idx = find(t_full >= last_beat_start);
    
    % Resetujemy wektor czasu, żeby ostatnie uderzenie zaczynało się znów od ~0 ms
    t_last = t_full(idx) - last_beat_start;
    V_last = V_full(idx);
    
    % Zapewnienie, że zwracany wektor pokrywa się DOKŁADNIE z okienkiem STOPTIME 
    % (wymagane, by funkcje normalizujące w Twoim UI nie zwariowały)
    idx_stop = find(t_last <= STOPTIME);
    
    t = t_last(idx_stop);
    V = V_last(idx_stop);

end