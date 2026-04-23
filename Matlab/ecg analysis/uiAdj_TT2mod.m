clearvars

% Parametry pacjenta
patient = 'normal_young_male'; offset = 169;

% --- ODCZYT CT (Potrzebne dla wrappera) ---
[~, ~, ventri_epi_ver_idx] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle_epi.tri'));
[~, ~, ventri_endo_ver_idx] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle_endo.tri'));

CT = ones(length(ventri_epi_ver_idx) + length(ventri_endo_ver_idx), 1);
for i = 1:length(ventri_endo_ver_idx)
    CT(ventri_endo_ver_idx(i)) = 3;
end
% -----------------------------------------

% Konfiguracja suwaków: Nazwa, Zakres [min max], Wartość domyślna
sConfigs(1) = struct('Name', 'Phase 1', 'Limits', [0 2], 'Default', 1.0);
sConfigs(2) = struct('Name', 'Phase 2', 'Limits', [0.5 1.5], 'Default', 1.0);
sConfigs(3) = struct('Name', 'Phase 3', 'Limits', [0.5 1.5], 'Default', 1.0);

% Uruchomienie nowego UI. Używamy funkcji anonimowej, by podać CT do wrappera!
ui = TT2modUiAdjuster(@(pE, pN, d, r, L, A) tt2_wrapper(pE, pN, d, r, L, A, CT), patient, offset, [], sConfigs);

% Funkcja opakowująca model TT2mod dla nowej wersji UI
function [tmp, sig_sim] = tt2_wrapper(pEpi, pEndo, dep, rep, L, A, CT)
    HT = 0.1; STOPTIME = 500; threshold = -40;
    
    % Generowanie dwóch szablonów AP
    [tE, VE] = TenTusscher2mod(HT, STOPTIME, 1, pEpi);
    [tN, VN] = TenTusscher2mod(HT, STOPTIME, 3, pEndo);
    
    % Wyznaczanie punktów charakterystycznych (dep/rep) dla szablonów
    find_fid = @(t,V) [t(find(V>=threshold,1)), t(find(V(find(V>=threshold,1):end)<=threshold,1)+find(V>=threshold,1)-1)];
    fidE = find_fid(tE, VE);
    fidN = find_fid(tN, VN);
    
    % Normalizacja szablonów do 0-1
    VE = (VE - min(VE)) / (max(VE) - min(VE));
    VN = (VN - min(VN)) / (max(VN) - min(VN));
    
    N = length(dep);
    tmp = zeros(N, L);
    T_target = 1:L;
    
    % Upewnijmy się, że dep i rep są wektorami kolumnowymi do obliczeń macierzowych
    dep = dep(:);
    rep = rep(:);
    
    % --- WEKTORYZACJA EPIKARDIUM ---
    idxE = find(CT == 1);
    if ~isempty(idxE)
        scaleE = (fidE(2) - fidE(1)) ./ (rep(idxE) - dep(idxE));
        shiftE = fidE(1) - scaleE .* dep(idxE);
        
        % Macierz czasów (rozmiar: length(idxE) x L)
        T_mappedE = scaleE .* T_target + shiftE;
        
        % Przycięcie czasów poza zakresem, by uniknąć błędu ekstrapolacji
        T_mappedE(T_mappedE < tE(1)) = tE(1);
        T_mappedE(T_mappedE > tE(end)) = tE(end);
        
        tmp(idxE, :) = interp1(tE, VE, T_mappedE, 'linear');
    end
    
    % --- WEKTORYZACJA ENDOKARDIUM ---
    idxN = find(CT == 3);
    if ~isempty(idxN)
        scaleN = (fidN(2) - fidN(1)) ./ (rep(idxN) - dep(idxN));
        shiftN = fidN(1) - scaleN .* dep(idxN);
        
        % Macierz czasów (rozmiar: length(idxN) x L)
        T_mappedN = scaleN .* T_target + shiftN;
        
        % Przycięcie czasów poza zakresem
        T_mappedN(T_mappedN < tN(1)) = tN(1);
        T_mappedN(T_mappedN > tN(end)) = tN(end);
        
        tmp(idxN, :) = interp1(tN, VN, T_mappedN, 'linear');
    end
    
    % Mnożenie przez macierz przejścia (lead field)
    sig_sim = A * tmp;
end