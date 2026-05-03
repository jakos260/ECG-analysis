patient = 'normal_young_male'; offset = 169; % cut p-wave
% patient = 'normal_male'; offset = 1;

dep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.dep'));
rep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.rep'));

A       = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2standard_12.mat'));
sig_ref = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/standard_12.refECG'));
Y = sig_ref(:, offset:end);

[~, ~, ventri_epi_ver_idx] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle_epi.tri'));
[~, ~, ventri_endo_ver_idx] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle_endo.tri'));

CT = ones(length(ventri_epi_ver_idx) + length(ventri_endo_ver_idx), 1);
for i = 1:length(ventri_endo_ver_idx)
    CT(ventri_endo_ver_idx(i)) = 3;
end
% ZAKŁADAMY, ŻE POSIADASZ W PRZESTRZENI ROBOCZEJ:
% Y - Macierz EKG (12 x 532)
% A - Macierz przekształceń (12 x 1500)
% cell_types - Wektor (1500 x 1) zawierający typ komórki dla każdego węzła (1=epi, 3=endo)

%% KROK 1: REKONSTRUKCJA ELEKTROGRAMÓW (Zewnątrzkomórkowych)
% Używamy regularyzacji Tichonowa 0-go rzędu, aby z EKG uzyskać sygnały na sercu.

lambda = 1e-2; % Parametr regularyzacji (do tuningu - zacznij od małych wartości)
num_nodes = size(A, 2); % 1500
time_samples = size(Y, 2); % 532

% Odwrócenie macierzy:
L = eye(num_nodes);
X_egm = (A' * A + (lambda^2) * (L' * L)) \ (A' * Y); 
% X_egm ma wymiar 1500 x 532 i zawiera wirtualne elektrogramy

%% KROK 2: EKSTRAKCJA CZASÓW AT i RT (Feature Extraction)
% Dla każdego węzła szukamy markerów fizjologicznych na pochodnej sygnału.

DT = zeros(num_nodes, 1);
RT = zeros(num_nodes, 1);
APD = zeros(num_nodes, 1); % Action Potential Duration

% Zdefiniuj okna czasowe (w próbkach), aby algorytm nie szukał RT w zespole QRS.
% Zakładam próbkowanie np. 1 próba = 1 ms.
window_QRS = 1:110;     % Okno szukania depolaryzacji
window_Twave = 130:500; % Okno szukania repolaryzacji

for i = 1:num_nodes
    egm_signal = X_egm(i, :);
    dVdt = diff(egm_signal); % Pierwsza pochodna przestrzenna sygnału
    
    % AT: Czas największego, najszybszego spadku potencjału w oknie QRS
    [~, at_idx] = min(dVdt(window_QRS));
    DT(i) = window_QRS(at_idx);
    
    % RT: Czas najszybszego wzrostu w oknie załamka T
    [~, rt_idx] = max(dVdt(window_Twave));
    RT(i) = window_Twave(rt_idx);
    
    % Wyliczenie wymaganego czasu trwania potencjału dla danego węzła
    APD(i) = RT(i) - DT(i); 
end

%% KROK 3: GENERACJA MAPY POTENCJAŁÓW CZYNNOŚCIOWYCH (AP)
% Używamy Twojego modelu Ten Tusscher'a, uwzględniając czasy aktywacji (AT).

AP_matrix = zeros(num_nodes, time_samples) - 85; % Inicjalizacja potencjałem spoczynkowym (-85 mV)
HT = 1; % Krok całkowania (zakładam 1 ms dla uproszczenia z macierzą Y)
phases_mod = [1, 1, 1]; % Startowe modyfikatory prądów faz 1, 2, 3

for i = 1:num_nodes    
    % Twój generator zakłada start AP w czasie t=0. 
    % Będziemy musieli wygenerować sygnał i przesunąć go w czasie.
    STOPTIME_AP = time_samples - AT(i) + 1; % Ile czasu zostało do końca symulacji od momentu AT
    
    if STOPTIME_AP > 0
        % Generacja "pudełka" potencjału czynnościowego
        [t, V] = TenTusscher2mod(HT, STOPTIME_AP, CT(i), phases_mod);
        
        % Mapowanie wygenerowanego sygnału w odpowiednie miejsce macierzy wynikowej
        % Wypełniamy od momentu aktywacji (AT) do końca okna
        length_to_insert = min(length(V), STOPTIME_AP);
        AP_matrix(i, AT(i):(AT(i)+length_to_insert-1)) = V(1:length_to_insert);
    end
end

% Na koniec AP_matrix (1500x532) posłuży jako wejście do macierzy transformacji A, 
% aby wygenerować EKG_sim i porównać z prawdziwym EKG!

%%
ECG_sim = A * AP_matrix;
leadv12( ...
    Y(:,1:500), ...
    ECG_sim(:,1:500), ...
    'paperspeed', 200, ...
    'amplification', 1 ...
    )