function LUT = getTmpLut_ResizedApd(start_apd, stop_apd, step_apd)
    % GETTMPLUT_RESIZEDAPD Generuje słownik potencjałów czynnościowych (LUT)
    % poprzez wygenerowanie jednego bazowego AP i jego czasową interpolację
    % do pożądanych szerokości (50% - 50% amplitudy).
    
    disp('Generowanie bazowego potencjału czynnościowego (AP)...');
    
    % 1. Parametry symulacji bazowej
    dt = 0.1; % 1 ms = 10 próbek (próbkowanie 10 kHz)
    sim_time = 1500; % Wydłużony czas (1500 ms), aby przy rozciąganiu AP nie obciąć sygnału
    
    % Wywołanie modelu bazowego dla sztywnego modyfikatora fazy 3: [1, 1, 1]
    [t_base, V_base] = wrapper_TenTusscher2mod(dt, sim_time, 1, [1, 1, 1], 150);
    
    % Normalizacja do zakresu [0, 1]
    V_norm = (V_base - min(V_base)) / (max(V_base) - min(V_base));
    
    % 2. Wyznaczenie bazowego APD (od 50% narastania do 50% opadania)
    idx_dep = find(V_norm >= 0.5, 1, 'first');
    idx_rep = find(V_norm(idx_dep:end) <= 0.5, 1, 'first') + idx_dep - 1;
    
    if isempty(idx_rep)
        error('Błąd: Bazowy AP nie zrepolaryzował się do 50% w zadanym czasie symulacji.');
    end
    
    t_dep_base = t_base(idx_dep);
    base_apd = t_base(idx_rep) - t_dep_base;
    
    disp(['Bazowy AP wygenerowany. Rzeczywiste APD (50%-50%): ', num2str(base_apd), ' ms.']);
    
    % 3. Przygotowanie wektora docelowych APD
    step_val = abs(step_apd);
    if start_apd <= stop_apd
        target_apds = start_apd : step_val : stop_apd;
    else
        target_apds = start_apd : -step_val : stop_apd;
    end
    
    num_targets = length(target_apds);
    
    % Prealokacja struktury dla wydajności
    LUT = repmat(struct('phase_mod', 1, 'V', [], 't', [], 't_dep', [], 'raw_APD', [], 'APD', []), 1, num_targets);
    
    disp('Skalowanie czasowe bazowego AP do docelowych szerokości...');
    
    % 4. Skalowanie i tworzenie finalnego LUT
    for k = 1:num_targets
        target = target_apds(k);
        
        % Współczynnik skalowania (proporcja docelowego APD do bazowego)
        scale_factor = target / base_apd;
        
        % Transformacja osi czasu. Zakotwiczamy punkt t_dep_base, 
        % aby nie przesunąć momentu aktywacji w czasie (upstroke pozostaje na miejscu).
        t_stretched = t_dep_base + (t_base - t_dep_base) * scale_factor;
        
        % Resampling (interpolacja) sygnału do oryginalnej, stałej siatki czasu (10 kHz)
        % Używamy wartości 'V_norm(end)' do ekstrapolacji po prawej stronie (aby uniknąć NaN)
        V_resized = interp1(t_stretched, V_norm, t_base, 'linear', V_norm(end));
        
        % Korekta błędów numerycznych z interpolacji
        V_resized(V_resized < 0) = 0;
        V_resized(V_resized > 1) = 1;
        
        % Odszukanie nowych indeksów do weryfikacji ostatecznej szerokości
        idx_dep_new = find(V_resized >= 0.5, 1, 'first');
        idx_rep_new = find(V_resized(idx_dep_new:end) <= 0.5, 1, 'first') + idx_dep_new - 1;
        
        if ~isempty(idx_rep_new)
            raw_apd = t_base(idx_rep_new) - t_base(idx_dep_new);
        else
            raw_apd = NaN;
        end
        
        % Zapisanie wyników
        LUT(k).phase_mod = 1; 
        LUT(k).V = V_resized;
        LUT(k).t = t_base;
        LUT(k).t_dep = t_base(idx_dep_new);
        LUT(k).raw_APD = raw_apd; 
        LUT(k).APD = target;
        
        % Opcjonalny pasek postępu, jeśli wykorzystujesz osobną funkcję w swoim workspace
        % printProgress(k, num_targets, sprintf('Tworzenie APD: %.1f ms', target));
    end
    
    disp(['Zakończono. Utworzono ', num2str(num_targets), ' idealnie dopasowanych szablonów.']);
end