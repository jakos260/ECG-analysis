function [t, V] = wrapper_TenTusscher2mod(HT, STOPTIME, CT, phases_mod, HR, ISO)
% WRAPPER_TENTUSSCHER2MOD Wrapper dla modelu Ten Tusscher 2
% Zwraca dokładnie 1 przebieg potencjału czynnościowego.

    % --- Obsługa opcjonalnych argumentów ---
    % Jeśli przekazano mniej niż 5 argumentów lub HR jest puste, ustaw domyślnie 60
    if nargin < 5 || isempty(HR)
        HR = 60;
    end

    % Jeśli przekazano mniej niż 6 argumentów lub ISO jest puste, ustaw domyślnie 0
    if nargin < 6 || isempty(ISO)
        ISO = 0;
    end

    % --- Parametry gwarantujące dokładnie JEDEN potencjał (1 AP) ---
    Stim_I = -52;   % Amplituda impulsu [A/F] (skalar = 1 impuls)
    Stim_T = 1;     % Czas trwania impulsu [ms]
    
    % Ustawiamy Stim_Int na niewielką wartość (np. 1 ms). 
    % Dzięki temu pojedynczy impuls pojawi się niemal od razu na początku 
    % okna czasowego (t=1), a symulacja potrwa aż do osiągnięcia STOPTIME.
    Stim_Int = 1;   

    % --- Wywołanie głównej funkcji modelu ---
    [t, V] = TenTusscher2mod(HT, STOPTIME, Stim_I, Stim_T, ISO, Stim_Int, CT, phases_mod);

end