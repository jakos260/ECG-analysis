%% ŁADOWANIE DANYCH 
patient = 'normal_young_male'; 
offset = 169; 
STOPTIME = 500; % Przykładowa długość sygnału

% Zamiast invInit, tworzymy strukturę ręcznie:
GEOM = struct();
GEOM.subject = patient;
GEOM.beat = 'beat1';
GEOM.type = 'ventricles';

% Prawdziwe czasy (Ground Truth) - idealne do późniejszego porównania z wynikiem inverse!
true_dep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.dep'));
true_rep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.rep'));

% Geometria serca
[GEOM.VER, GEOM.ITRI] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.tri'));

% Macierz transferu (Ventricles -> 12-lead ECG)
GEOM.AMA = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2BSM_(nijmegen_64).mat'));

% Sygnały EKG (Twoje "BSM" to teraz po prostu 12 odprowadzeń)
sig_ref = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/BSM_(nijmegen_64).refECG'));
GEOM.BSM = sig_ref(:, offset:offset+STOPTIME-1);
GEOM.TVER = zeros(size(GEOM.BSM, 1), 3);

% Podstawowa macierz odległości (GEOM.DIST)
GEOM.DIST = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.dst3d'));

% Macierz odległości po powierzchni (GEOM.DISTsurf)
GEOM.DISTsurf = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.dst2d'));

% Macierze z wagami anizotropowymi (GEOM.DIST2W oraz GEOM.ADJ2W)
GEOM.DIST2W = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.dstanis'));
GEOM.ADJ2W = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.adjanis'));

% Liczba węzłów na podstawie wczytanej wcześniej siatki ventricle.tri
num_nodes = length(GEOM.VER);

GEOM.purkinjever = zeros(num_nodes, 1);
GEOM.Lpurkinjever = zeros(num_nodes, 1);
GEOM.Rpurkinjever = zeros(num_nodes, 1);
GEOM.Rfreewallver = zeros(num_nodes, 1);
GEOM.endoVER = zeros(num_nodes, 1);

GEOM.neigh = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.adjneigh'));

% Pusta zmienna pS, która jest podawana w funkcji getSmode
GEOM.pS = [];
GEOM.ADJ = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.adj3d'));
GEOM.anisotropyRatio = 0.5; % value only to save in logs

% GEOM.specs(2) - początek analizowanego sygnału (np. początek QRS)
% GEOM.specs(3) - koniec QRS (na tej podstawie liczone jest qrsduration)
% GEOM.specs(4) - używane do estymacji repolaryzacji
% GEOM.specs(5) - koniec analizowanego sygnału (koniec załamka T)
GEOM.specs = [0, 1, 100, 200, STOPTIME];

% stara struktura GEOM.SPECS
GEOM.SPECS = struct();
GEOM.SPECS.onsetqrs = GEOM.specs(2);
GEOM.SPECS.onsetp = 0; % used only in atria mode
GEOM.SPECS.qrstduration = GEOM.specs(5) - GEOM.specs(2) + 1;
GEOM.SPECS.time_Jpoint = GEOM.specs(3);
GEOM.SPECS.time_Vstim = 1; % used only in stim mode
GEOM.SPECS.time_apexT = 275;
GEOM.SPECS.depSlope = 1.0;
GEOM.SPECS.repCorrection = 0;
GEOM.SPECS.endtwave = 400;
GEOM.SPECS.useCumsum = false;
GEOM.SPECS.initialSlope = 1;
GEOM.SPECS.plateauslope = 0.02;
GEOM.SPECS.repslope = 0.045;
GEOM.SPECS.scaleAmpl = 15;

% Opcjonalnie: uśrednienie sygnału (baseline correction)
% GEOM.BSM = zeromean(GEOM.BSM);

GEOM.LAY = lay_nim64();

%% BUDOWA SŁOWNIKA LUT (Ten Tusscher 2) I ŁADOWANIE KLASTRÓW
disp('Ładowanie podziału na segmenty (AHA)...');
GEOM.RegionIdx = load(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles_labels_of_20.segments'));

disp('Budowanie słownika LUT dla fizjologicznych kształtów APD...');
phases_sweep = 0.5:0.1:1.5;
LUT = struct();

for i = 1:length(phases_sweep)
    % Uruchomienie Twojego generatora dla danej modyfikacji fazy 3
    [t_tt2, V_tt2] = wrapper_TenTusscher2mod(0.1, 500, 1, [1, 1, phases_sweep(i)], 60, 0);
    
    % Normalizacja potencjału dla macierzy AMA
    V_norm = (V_tt2 - min(V_tt2)) / (max(V_tt2) - min(V_tt2));
    
    % Wyliczenie dokładnych markerów czasowych dla tego wzorca
    idx_dep = find(V_norm >= 0.5, 1, 'first');
    idx_rep = find(V_norm(idx_dep:end) <= 0.1, 1, 'first') + idx_dep - 1;
    if isempty(idx_rep), idx_rep = length(V_norm); end
    
    LUT(i).phase_mod = phases_sweep(i);
    LUT(i).V = V_norm;
    LUT(i).t = t_tt2;
    LUT(i).t_dep = t_tt2(idx_dep);
    LUT(i).APD = t_tt2(idx_rep) - LUT(i).t_dep;
end
GEOM.LUT = LUT; % Zapisujemy do struktury przekazywanej do invfunc

%% ESTYMACJA POCZĄTKOWA I PROBLEM ODWROTNY
initialvelocity = 0.4;
mudep = 5e-5;
murep = 1e-4;

% 1. Szukanie ognisk początkowych (Initial Estimate)
disp('Uruchamiam multifociscan...');
[measinit.foci, measinit.dep, measinit.outp] = multifociscan(GEOM, 1, 0);

%%
measinit.rep = initRep(GEOM, measinit.dep);

% 2. Dokładna optymalizacja (Inverse Problem)
disp('Uruchamiam inversefunc...');
meas = my_inversefunc(GEOM, measinit.dep, measinit.rep, mudep, murep);

%% 4. PORÓWNANIE Z ORYGINAŁEM
error_dep = abs(meas.depfinal - true_dep);
error_rep = abs(meas.repfinal - true_rep);
disp(['Średni błąd depolaryzacji: ', num2str(mean(error_dep)), ' ms']);
disp(['Średni błąd repolaryzacji: ', num2str(mean(error_rep)), ' ms']);


%% 5. Qtriplot
q = initQtripy();

q.reset();
q.disable_debounce();
q.set_panels_number(2, 1);

depolarisation time error
q.set_active_panel(2, 1);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(0.3);
q.values(error_dep);
q.gradient_bins(10);
q.text("error_dep", [0.6, 0.85]);

simulated repolarisation time
q.set_active_panel(1, 1);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(0.3);
q.values(error_rep);
q.gradient_bins(10);
q.text("error_rep", [0.1, 0.85]);