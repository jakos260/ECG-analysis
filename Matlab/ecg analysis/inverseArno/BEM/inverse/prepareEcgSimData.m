%% 1. ŁADOWANIE DANYCH 
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
GEOM.AMA = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2standard_12.mat'));

% Sygnały EKG (Twoje "BSM" to teraz po prostu 12 odprowadzeń)
sig_ref = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/standard_12.refECG'));
GEOM.BSM = sig_ref(:, offset:offset+STOPTIME-1);

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

% Pusta zmienna pS, która jest podawana w funkcji getSmode
GEOM.pS = [];

%% 2. RĘCZNE DEFINIOWANIE MARKERÓW (zamiast prepareECG)
% GEOM.specs(2) - początek analizowanego sygnału (np. początek QRS)
% GEOM.specs(3) - koniec QRS (na tej podstawie liczone jest qrsduration)
% GEOM.specs(4) - używane do estymacji repolaryzacji
% GEOM.specs(5) - koniec analizowanego sygnału (koniec załamka T)
GEOM.specs = [0, 1, 100, 200, STOPTIME];

% Opcjonalnie: uśrednienie sygnału (baseline correction)
% GEOM.BSM = zeromean(GEOM.BSM);

%% 3. ESTYMACJA POCZĄTKOWA I PROBLEM ODWROTNY
initialvelocity = 0.4;
mudep = 5e-6;
murep = 1e-5;

% 1. Szukanie ognisk początkowych (Initial Estimate)
disp('Uruchamiam multifociscan...');
[measinit.foci, measinit.dep, measinit.outp] = multifociscan(GEOM, 1, 0);

measinit.rep = initRep(GEOM, measinit.dep);

% 2. Dokładna optymalizacja (Inverse Problem)
disp('Uruchamiam inversefunc...');
meas = inversefunc(GEOM, measinit.dep, measinit.rep, ...
    'repopt', 'apd', 'maxiter', 400, 'mudep', mudep, 'murep', murep);

%% 4. PORÓWNANIE Z ORYGINAŁEM
dep_error = abs(meas.depfinal - true_dep);
disp(['Średni błąd depolaryzacji: ', num2str(mean(dep_error)), ' ms']);