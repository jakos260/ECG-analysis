%%
I0  = -24;  % Stimulus strength (A/F)
T0  = 1;    % Pulse duration (ms)
RR0 = 1000; % Standard RR interval (ms)

% Create pulse train

% Normalized stimulus current for each pulse:
I = [1 1 1 1 1 1 1 1 1 1 1 1 1 ];
% Normalized pulse duration for each pulse:
T = [1 1 1 1 1 1 1 1 1 1 1 1 1 ];
% Normalized interpulse interval after start of each pulse
RR = [0.1 1 1 1 1 1 1 1 1 1 1 1 1 ];
% ( RR(1) is the starttime of the first pulse)

% modify phases
phases_mod = [1, 1, 1];
modifications = [0, 0.8, 1, 2];

cell_types = [1, 3]; % 1 - Epicardial, 3 - Endocardial
cell_names = {'Epicardial', '', 'Endocardial'};

for CellType = cell_types
    figure('Name', cell_names{CellType}, 'NumberTitle', 'off');
    sgtitle(['Model Ten Tusscher 2: ', cell_names{CellType}]); % figure title
    
    for phase_idx = 1:length(phases_mod)
        subplot(3, 1, phase_idx);
        hold on;
        
        for mod_idx = 1:length(modifications)
            
            % Resetowanie phases_mod do wartości bazowych [1, 1, 1]
            current_phases_mod = [1, 1, 1];
            % Podstawienie konkretnej modyfikacji pod odpowiednią fazę
            current_phases_mod(phase_idx) = modifications(mod_idx);
            
            % Wywołanie modelu
            [t, VAP] = TenTusscher2mod(I0*I, T0*T, RR0*RR, CellType, current_phases_mod);
            
            % Ekstrakcja sygnału
            n1 = length(t) - 55000;
            n2 = length(t) - 25000;
            t_singleAP = t(n1:n2) - t(n1);
            V_singleAP = VAP(n1:n2);
            
            % Rysowanie wykresu z dodaniem nazwy do legendy
            plot(t_singleAP, V_singleAP, 'LineWidth', 1.5, ...
                 'DisplayName', ['Mod = ', num2str(modifications(mod_idx))]);
        end
        
        % Oznaczenia subplota
        title(['Modification of phase ', num2str(phase_idx)]);
        xlabel('t (ms)');
        ylabel('V_{membrane} (mV)');
        legend('Location', 'best');
        grid on;
        hold off;
    end
end

%% 
% ZADANE ZMIENNE (przykład, podmień na swoje dane, jeśli trzeba)
% dep = rand(100, 1) * 20;       % Czasy depolaryzacji w ms
% rep = dep + 250 + rand(100, 1)*50; % Czasy repolaryzacji w ms
% isEpi = randi([0, 1], 100, 1); % 1 - Epi, 0 - Endo
N_nodes = length(dep);

% 1. PARAMETRY CZASOWE MACIERZY WYNIKOWEJ
dt = 1; % Krok czasowy dla macierzy wyjściowej w ms
t_max = max(rep) + 200; % Symulacja do najpóźniejszej repolaryzacji + bufor
t_global = 0:dt:t_max;
TMP_matrix = zeros(N_nodes, length(t_global)); % Prealokacja macierzy (Ważne dla wydajności!)

% 2. GENEROWANIE TEMPLATE'ÓW (BEZ MODYFIKACJI)
phases_mod = [1, 1, 1];

% --- Template Epicardial (CellType = 1) ---
[t_e, V_e] = TenTusscher2mod(I0*I, T0*T, RR0*RR, 1, phases_mod);
n1 = length(t_e)-55000; n2 = length(t_e)-25000;
t_tmpl_epi = t_e(n1:n2) - t_e(n1);
V_tmpl_epi = V_e(n1:n2);

% --- Template Endocardial (CellType = 3) ---
[t_en, V_en] = TenTusscher2mod(I0*I, T0*T, RR0*RR, 3, phases_mod);
t_tmpl_endo = t_en(n1:n2) - t_en(n1);
V_tmpl_endo = V_en(n1:n2);

% 3. ZNALEZIENIE PUNKTÓW REFERENCYJNYCH (-40 mV) W TEMPLATE'ACH
threshold = -40;

% Funkcja pomocnicza do szukania indeksów depolaryzacji i repolaryzacji
find_fiducials = @(t, V) deal(...
    t(find(V >= threshold, 1, 'first')), ... % Czas depolaryzacji
    t(find(V(find(V >= threshold, 1, 'first'):end) <= threshold, 1, 'first') + find(V >= threshold, 1, 'first') - 1) ... % Czas repolaryzacji
);

[t_dep_epi, t_rep_epi] = find_fiducials(t_tmpl_epi, V_tmpl_epi);
[t_dep_endo, t_rep_endo] = find_fiducials(t_tmpl_endo, V_tmpl_endo);

% 4. WYPEŁNIANIE MACIERZY (Time-Warping dla każdego węzła)
disp('Rozpoczęto mapowanie potencjałów do węzłów...');
for i = 1:N_nodes
    % Wybór odpowiedniego template'u na podstawie isEpi
    if isEpi(i) == 1
        t_tmpl = t_tmpl_epi; V_tmpl = V_tmpl_epi;
        t_dep_T = t_dep_epi; t_rep_T = t_rep_epi;
    else
        t_tmpl = t_tmpl_endo; V_tmpl = V_tmpl_endo;
        t_dep_T = t_dep_endo; t_rep_T = t_rep_endo;
    end
    
    dep_i = dep(i);
    rep_i = rep(i);
    
    % Wektor zmapowanego czasu dla interpolacji
    t_mapped = zeros(size(t_global));
    
    % Faza 1: Przed depolaryzacją (zachowujemy oryginalne tempo)
    idx_before = t_global < dep_i;
    t_mapped(idx_before) = t_global(idx_before) - dep_i + t_dep_T;
    
    % Faza 2: Czas trwania potencjału czynnościowego (AP) - SKALOWANIE
    idx_ap = (t_global >= dep_i) & (t_global <= rep_i);
    % Współczynnik skalowania: oryginalny APD / docelowy APD
    scale = (t_rep_T - t_dep_T) / (rep_i - dep_i); 
    t_mapped(idx_ap) = t_dep_T + (t_global(idx_ap) - dep_i) * scale;
    
    % Faza 3: Po repolaryzacji (zachowujemy oryginalne tempo)
    idx_after = t_global > rep_i;
    t_mapped(idx_after) = t_rep_T + (t_global(idx_after) - rep_i);
    
    % Ograniczenie brzegowe, żeby uniknąć wyjścia poza zakres t_tmpl (utrzymanie potencjału spoczynkowego)
    t_mapped = max(min(t_mapped, t_tmpl(end)), t_tmpl(1));
    
    % Interpolacja kształtu - używamy 'pchip' aby uniknąć przesterowań (overshoots)
    TMP_matrix(i, :) = interp1(t_tmpl, V_tmpl, t_mapped, 'pchip');
end
disp('Mapowanie zakończone!');

% --- (Opcjonalnie) Wyrysowanie 3 losowych węzłów dla weryfikacji ---
figure;
hold on;
rand_nodes = randi(N_nodes, 3, 1);
colors = lines(3);
for k = 1:3
    node = rand_nodes(k);
    plot(t_global, TMP_matrix(node, :), 'Color', colors(k,:), ...
        'DisplayName', sprintf('Węzeł %d (Epi: %d, Dep: %.1f, Rep: %.1f)', node, isEpi(node), dep(node), rep(node)));
    % Zaznaczenie punktów depolaryzacji i repolaryzacji
    plot(dep(node), -40, 'o', 'Color', colors(k,:), 'HandleVisibility', 'off');
    plot(rep(node), -40, 'x', 'Color', colors(k,:), 'HandleVisibility', 'off');
end
yline(-40, 'k--', 'Próg -40mV', 'HandleVisibility', 'off');
xlabel('Czas (ms)');
ylabel('V_{membrane} (mV)');
title('Zmapowane przebiegi TMP dla wybranych węzłów');
legend('Location', 'best');
grid on;
hold off;