%% ŁADOWANIE DANYCH median
subject = 'normal_young_male'; 
offset = 169; 
STOPTIME = 500; % Przykładowa długość sygnału

% Zamiast invInit, tworzymy strukturę ręcznie:
GEOM = struct();
GEOM.subject = subject;
GEOM.beat = 'beat1';
GEOM.type = 'ventricles';

% Prawdziwe czasy (Ground Truth) - idealne do późniejszego porównania z wynikiem inverse!
true_dep = loadmat(append(DATA_PATH, 'other/ECGsim_data/', subject, '/ventricular_beats/beat1/user.dep'));
true_rep = loadmat(append(DATA_PATH, 'other/ECGsim_data/', subject, '/ventricular_beats/beat1/user.rep'));

% Geometria serca
[GEOM.VER, GEOM.ITRI] = loadtri_ecgsim(append(DATA_PATH, 'other/ECGsim_data/', subject, '/model/ventricle.tri'));

% Macierz transferu (Ventricles -> 12-lead ECG)
GEOM.AMA = loadmat(append(DATA_PATH, 'other/ECGsim_data/', subject, '/model/ventricles2BSM_(nijmegen_64).mat'));

% Sygnały EKG (Twoje "BSM" to teraz po prostu 12 odprowadzeń)
sig_ref = loadmat(append(DATA_PATH, 'other/ECGsim_data/', subject, '/ecgs/BSM_(nijmegen_64).refECG'));
GEOM.BSM = sig_ref(:, offset:offset+STOPTIME-1);
GEOM.TVER = zeros(size(GEOM.BSM, 1), 3);

% Podstawowa macierz odległości (GEOM.DIST)
GEOM.DIST = loadmat(append(DATA_PATH, 'other/ECGsim_data/', subject, '/model/ventricle.dst3d'));

% Macierz odległości po powierzchni (GEOM.DISTsurf)
GEOM.DISTsurf = loadmat(append(DATA_PATH, 'other/ECGsim_data/', subject, '/model/ventricle.dst2d'));

% Macierze z wagami anizotropowymi (GEOM.DIST2W oraz GEOM.ADJ2W)
GEOM.DIST2W = loadmat(append(DATA_PATH, 'other/ECGsim_data/', subject, '/model/ventricle.dstanis'));
GEOM.ADJ2W = loadmat(append(DATA_PATH, 'other/ECGsim_data/', subject, '/model/ventricle.adjanis'));

% Liczba węzłów na podstawie wczytanej wcześniej siatki ventricle.tri
num_nodes = length(GEOM.VER);

GEOM.purkinjever = zeros(num_nodes, 1);
GEOM.Lpurkinjever = zeros(num_nodes, 1);
GEOM.Rpurkinjever = zeros(num_nodes, 1);
GEOM.Rfreewallver = zeros(num_nodes, 1);
GEOM.endoVER = zeros(num_nodes, 1);

GEOM.neigh = loadmat(append(DATA_PATH, 'other/ECGsim_data/', subject, '/model/ventricle.adjneigh'));

% Pusta zmienna pS, która jest podawana w funkcji getSmode
GEOM.pS = [];
GEOM.ADJ = loadmat(append(DATA_PATH, 'other/ECGsim_data/', subject, '/model/ventricle.adj3d'));
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
disp('Loading AHA segment definitions...');
GEOM.RegionIdx = load(append(DATA_PATH, 'other/ECGsim_data/', subject, '/model/ventricles.typ'));

% Tworzymy słownik z idealnie równymi krokami od 200 do 350 ms z krokiem co 1 ms
GEOM.LUT = getTmpLut_niceApd(200, 460, 1);

    q
%%
% (Initial Estimate - foci skan)
disp('Uruchamiam multifociscan...');
[measinit.foci, measinit.dep, measinit.outp] = multifociscan(GEOM, 1, 0);

%%
% initialvelocity = 0.4;
mudep = 1e-6;
murep = 2e-6;

measinit.rep = initRep(GEOM, measinit.dep);

disp('Uruchamiam inversefunc...');
meas = my_inversefunc(GEOM, measinit.dep, measinit.rep, mudep, murep);

%% old aproach
mudep = 5e-5;
murep = 1e-6;
estamp = 0;
maxit = 400;
mrd   = 0.25;
scanmode = 4;
wrd = 0;
blmode = 1;
lpass = 1;

meas_old = inversefunc( ...
    GEOM,measinit.dep,measinit.rep, ...
    'estimateampl',estamp, ...
    'repopt','apd', ...
    'maxiter',maxit, ...
    'mudep',mudep, ...
    'murep',murep, ...
    'minrd',mrd, ...
    'mode',scanmode,...
    'weighedrd',wrd, ...
    'blmode',blmode, ...
    'lpass',lpass ...
    );

%% 4. PORÓWNANIE Z ORYGINAŁEM
error_dep = abs(meas.depfinal - true_dep);
error_rep = abs(meas.repfinal - true_rep);
error_dep_old = abs(meas_old.depfinal - true_dep);
error_rep_old = abs(meas_old.repfinal - true_rep);

fprintf('Mean error dep=%.1f[ms] rep=%.1f[ms]\n', mean(error_dep), mean(error_rep));
fprintf('Mean error dep=%.1f[ms] rep=%.1f[ms] OLD\n', mean(error_dep_old), mean(error_rep_old));

% rel_diff_rep = RelativeDistance()

% data            = [];
% data.T          = T;
% data.meas       = meas;
% data.measinit   = measinit;
% plotecg(GEOM, data, S.scanmode, S.lpass)

%% 5. Qtriplot
transparency = 0.4;
gradient_bins = 6;

q = initQtripy();

q.reset();
q.disable_debounce();
q.set_panels_number(2, 2);

% depolarisation time error
q.set_active_panel(1, 1);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(error_dep);
q.gradient_bins(gradient_bins);
q.text(sprintf("error_dep av=%.0f[ms]", mean(error_dep)), [0.01, 0.95]);

% simulated repolarisation time
q.set_active_panel(2, 1);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(error_rep);
q.gradient_bins(gradient_bins);
q.text(sprintf("error_rep av=%.0f[ms]", mean(error_rep)), [0.51, 0.95]);


q.set_active_panel(1, 2);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(error_dep_old);
q.gradient_bins(gradient_bins);
q.text(sprintf("error_dep_old av=%.0f[ms]", mean(error_dep_old)), [0.01, 0.45]);

% simulated repolarisation time
q.set_active_panel(2, 2);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(error_rep_old);
q.gradient_bins(gradient_bins);
q.text(sprintf("error_rep_old av=%.0f[ms]", mean(error_rep_old)), [0.51, 0.45]);

%% ECGI Summit 2026

transparency = 0.4;
gradient_bins = 6;

q = initQtripy();

q.reset();
q.disable_debounce();
q.set_panels_number(1, 2);
q.background_color("white");

% depolarisation time error
q.set_active_panel(1, 1);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(meas.repfinal);
q.gradient_bins(gradient_bins);
q.text(sprintf("av=%.0f[ms]", mean(meas.repfinal)), [0.2, 0.95]);

% simulated repolarisation time
q.set_active_panel(1, 2);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(meas_old.repfinal);
q.gradient_bins(gradient_bins);
q.text(sprintf("av=%.0f[ms]", mean(meas_old.repfinal)), [0.2, 0.45]);


%%
N = 4;
% Smooth the signals (using a Gaussian filter with a small window size)
window_size = 5;

% Note: If offset:offset+400 contains 401 samples, ensure simulated_ecg indices match (e.g., 1:401).
ecg_ref = smoothdata(sig_ref(N, offset:offset+400), 'gaussian', window_size);
ecg_sim = smoothdata(meas.simulated_ecg(N, 1:400), 'gaussian', window_size); 
ecg_sim_old = smoothdata(meas_old.simulated_ecg(N, 1:400), 'gaussian', window_size);

figure(67);
clf; % Clear the figure to prevent overlapping when re-running
set(gcf, 'Color', 'w');

% First subplot: New simulation (LUT based)
subplot(2,1,1);
yline(0, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); % Add gray zero-line first so it stays in the background
hold on;
plot(ecg_ref, 'k', 'LineWidth', 2.0); % Thick black line for reference
plot(ecg_sim, 'r', 'LineWidth', 2.0); % Thick red line for simulation

% Format title to include the Relative Distance (RD) value
title(sprintf('Reference vs New LUT Optimization (RD = %.3f)', meas.rdfinal), 'FontSize', 12);
axis off; % Remove all axes, borders, and ticks

% Second subplot: Old simulation
subplot(2,1,2);
yline(0, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); % Add gray zero-line
hold on;
plot(ecg_ref, 'k', 'LineWidth', 2.0);
plot(ecg_sim_old, 'r', 'LineWidth', 2.0);

% Format title to include the Relative Distance (RD) value from the old model
title(sprintf('Reference vs Old Optimization (RD = %.3f)', meas_old.rdfinal), 'FontSize', 12);
axis off; % Remove all axes, borders, and ticks