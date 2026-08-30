% #############################################
% ### Computing in Cardiology 2026 - Madrid ###
% #############################################

num_subjects = 1; 

subject = sprintf('IKEM_Pat%03d', num_subjects);
path = 'C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data\Dataset\';
DATA = readGeomPeacsModelDataset(path, subject );

CineEcgResultPath = fullfile(path, subject, 'signals', sprintf('IKEM_Pat%03d.iecg', num_subjects));
EcgDataPath = fullfile(path, subject, 'signals', 'ECG_DATA');

N = 4; % 3 for LV pacing or 4 for RV pacing
cine_ecg = readCineEcgResults(CineEcgResultPath,'ecgdir', EcgDataPath, 'domedian'); 
BSM = cine_ecg.MEDIANDATA.VENTRICULAR{N}.beats.ECGbeat;

%%
STOPTIME = 500;

GEOM = struct();
GEOM.subject = DATA.subject;
GEOM.type = 'ventricles';

lead_system = 'x99Prague'; %'x65Nijmegen' 'x12plus3leads';
ref_signal = sprintf('s%d', N);

beat_no = 3;
beat_time =  DATA.VENTR.SIGNALS.INFO.(ref_signal).beats(beat_no);
beat_cut = [beat_time-200 beat_time+400];

GEOM.VER        = DATA.VENTR.geom.VER;
GEOM.ITRI       = DATA.VENTR.geom.ITRI;
GEOM.AMA        = DATA.VENTR.AMA.(lead_system);
GEOM.BSM        = BSM;
GEOM.DIST       = DATA.VENTR.DIST3D;
GEOM.DISTsurf   = DATA.VENTR.DISTsurf;
GEOM.DIST2W     = DATA.VENTR.DISTanis;
GEOM.ADJ2W      = DATA.VENTR.ADJanis;   % ventricles.adjanis; 
GEOM.neigh      = DATA.VENTR.ADJneigh;  % ventricle.adjneigh
GEOM.ADJ        = DATA.VENTR.ADJ3D;     % ventricle.adj3d
GEOM.LAY        = loadmat('C:\Users\Admin\Documents\Projects\ecg_project\Scripts\Matlab\ecg_analysis\inverseArno\BEM\inverse\mla\prague99.mla');
GEOM.RegionIdx  = DATA.GEOM.ventr.segments;

num_nodes           = length(GEOM.VER);
GEOM.TVER           = zeros(size(GEOM.BSM, 1), 3);
GEOM.purkinjever    = zeros(num_nodes, 1);
GEOM.Lpurkinjever   = zeros(num_nodes, 1);
GEOM.Rpurkinjever   = zeros(num_nodes, 1);
GEOM.Rfreewallver   = zeros(num_nodes, 1);
GEOM.endoVER        = zeros(num_nodes, 1);


GEOM.pS = [];
GEOM.anisotropyRatio = 0.5; % value only to save in logs

% GEOM.specs(2) - początek analizowanego sygnału (np. początek QRS)
% GEOM.specs(3) - koniec QRS (na tej podstawie liczone jest qrsduration)
% GEOM.specs(4) - używane do estymacji repolaryzacji
% GEOM.specs(5) - koniec analizowanego sygnału (koniec załamka T)
GEOM.specs = [0, 1, 100, 200, STOPTIME];
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
GEOM.LUT = getApLut_ResizedApd(150, 475, 1);


maxAmpl = round(max(max(abs(GEOM.BSM))));
figure(99);
clf;
sigplot(GEOM.BSM,'',GEOM.LAY,1.3/maxAmpl,'b',1,0); 
title(subject + " " + DATA.VENTR.SIGNALS.INFO.(ref_signal).name, 'Interpreter', 'none');


%% multifoci skan

velocity = 1.0; % m/s
init = simplyfocusscan(GEOM, velocity);
init_dep = init.dep;

% init = lineargradientscan(GEOM, [0, 1, -0.2], velocity); % dep vector
init = lineargradientscan(GEOM, [-1, 0, 0], 0.5);
init_rep = init.rep;
% init_dep = zeros(size(init.dep, 1), size(init.dep, 2)) + init.r_time;
% init_rep = zeros(size(init.rep, 1), size(init.rep, 2)) + init.t_time;

% segments_list = [1,2, 4, 5];
% for i=1:size(init_rep, 1)
%     % if ismember(DATA.GEOM.ventr.segments(i), segments_list)
%     % if DATA.GEOM.ventr.segments(i) > 17
%     %     init_rep(i) = init_rep(i) * 1.2;
%     % else
%     %     init_rep(i) = init_rep(i) * 0.8;
%     % end
% 
%     if DATA.GEOM.ventr.endoVER(i) == 0      % epi
%         init_rep(i) = init_rep(i) * 1.0;
%     elseif DATA.GEOM.ventr.endoVER(i) == 1  % endo RV
%         init_rep(i) = init_rep(i) * 1.2;
%     elseif DATA.GEOM.ventr.endoVER(i) == 2  % 2 endo LV
%         init_rep(i) = init_rep(i) * 1;     
%     end
% end

INV = struct();
INV.AMA   = GEOM.AMA;
INV.SPECS = GEOM.SPECS;
INV.RegionIdx = DATA.GEOM.ventr.segments;
INV.lpass = 1; % Default lowpass filter setting
INV.useWeighedRd = 0;
INV.LUT = GEOM.LUT;

INV.BSM = GEOM.BSM;
INV.PHIREF = INV.BSM;
INV.normphi = norm(INV.PHIREF, 'fro');

usetimes = size(INV.PHIREF, 2);
INV.T = ones(size(INV.AMA, 2), 1) * (0 : usetimes - 1);

[INV.REGOP, INV.REGOPREP] = calcREGOP(GEOM, 1);
INV.MuValues = [0, 0]; % Dummy mu values

% Pack initial times into a cell array and setup defaults for notch/ampl
params_cell = {init_dep, init_rep};
notchopt    = struct('pot', zeros(size(GEOM.VER,1),1));
amplopt     = struct('pot', ones(size(GEOM.VER,1),1));

% forward model
TST = gettres_v_nparams(INV, params_cell, notchopt, amplopt);


% RESULTS
plot_ecg_signals(INV.PHIREF, TST.PHIA, GEOM.LAY, TST.rd);


%%

% panel 2 - thorax różnica

q2 = initQtripy();
q2.reset();
q2.disable_debounce();
q2.background_color("white");

q2.text(sprintf("diff | t = %d ms", t_peak), [0.3, 0.95]);
q2.surface(thorax_ver, thorax_itri);
q2.transparency(0.0);
diff = pot_ref_full - pot_sim;
q2.values(diff);
q2.gradient_bins(10);
q2.minmax_values(-1,1);

q2.markers(DATA.VENTR.LEADPOS.x99Prague, 'black', 10);


%% 3D hearts
transparency = 0.0;
gradient_bins = 25;
q = initQtripy();

q.reset();
q.disable_debounce();
q.set_panels_number(1,2);
q.background_color("white");
q.color_range(0, max(init_rep));

% depolarisation time
q.set_active_panel(1, 1);
q.text(sprintf("dep av=%.0f[ms]", mean(init_dep)), [0.1, 0.95]);

q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(init_dep);
q.gradient_bins(gradient_bins);

% repolarisation time
q.set_active_panel(1, 2);
q.text(sprintf("rep av=%.0f[ms]", mean(init_rep)), [0.1, 0.45]);

q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(init_rep);
q.gradient_bins(gradient_bins);

% q.axis();

%% --- PORÓWNANIE ROZKŁADU POTENCJAŁÓW NA KLATCE (REPOLARYZACJA) ---

% 1. Identyfikacja momentu repolaryzacji (szczyt fali T)
% Zakładamy, że fala T znajduje się w drugiej połowie analizowanego sygnału
% 1. Identyfikacja momentu repolaryzacji
half_time = round(size(INV.PHIREF, 2) / 2);
rms_ref = rms(INV.PHIREF, 1);
[~, t_peak_offset] = max(rms_ref(half_time:end));
t_peak = half_time + t_peak_offset - 1;

% 2. Pobranie wektorów potencjałów dla znalezionego czasu
pot_ref_99 = INV.PHIREF(:, t_peak); % Potencjały referencyjne (99 odprowadzeń)
pot_sim = DATA.VENTR.THORAX * TST.S(:, t_peak); % Symulacja na pełnej klatce

% Wykonujemy błyskawiczną estymację potencjałów serca z 99 elektrod
% przy użyciu macierzy INV.AMA i regularyzacji Tichonowa
lambda_vis = 1e-2; 
pot_heart_est = (INV.AMA' * INV.AMA + lambda_vis * eye(size(INV.AMA, 2))) \ (INV.AMA' * pot_ref_99);

% Rzutujemy uzyskane potencjały serca na PEŁNĄ siatkę klatki piersiowej
pot_ref_full = DATA.VENTR.THORAX * pot_heart_est;
% ---------------------------------------

% Geometria klatki piersiowej
thorax_ver  = DATA.GEOM.thorax.VER;
thorax_itri = DATA.GEOM.thorax.ITRI;

% 3. Inicjalizacja qtripy z dwoma panelami
q3 = initQtripy();
q3.reset();
q3.disable_debounce();
q3.set_panels_number(3, 1); % Dwa panele obok siebie (X=2, Y=1)
q3.background_color("white");
q3.text(sprintf("t = %d ms", t_peak), [0.45, 0.15]);

% Wspólna skala kolorów
max_val = max(max(abs(pot_ref_full)), max(abs(pot_sim)));
q3.color_range(-max_val, max_val);

% 4. Panel 1 (Lewy): Referencja
q3.set_active_panel(1, 1);
% Zmiana: wyświetlamy t_peak jako czas, a nie średnią wartość napięcia
q3.text(sprintf("ref"), [0.1, 0.95]);
q3.surface(thorax_ver, thorax_itri);
q3.transparency(0.0);
q3.values(pot_ref_full); % TERAZ WYSYŁAMY PEŁNĄ KLATKĘ, A NIE 99 PUNKTÓW
q3.gradient_bins(25);

for i = 1:size(DATA.VENTR.LEADPOS.x99Prague, 1)
    q3.markers(DATA.VENTR.LEADPOS.x99Prague, 'black', 10);
end

% 5. Panel 2 (Prawy): Symulacja
q3.set_active_panel(2, 1);
q3.text(sprintf("sim"), [0.4, 0.95]);
q3.surface(thorax_ver, thorax_itri);
q3.transparency(0.0);
q3.values(pot_sim);
q3.gradient_bins(25);

% panel 3 - różnica
q3.set_active_panel(3, 1);
q3.text(sprintf("diff"), [0.7, 0.95]);
q3.surface(thorax_ver, thorax_itri);
q3.transparency(0.0);
q3.values(pot_ref_full - pot_sim);
q3.gradient_bins(25);


%% patterns


close all

% time setup
HT = 0.1;                       % timestep (ms)
STOPTIME = 600; % Duration of the simulation (ms)

% data read
patient = 'normal_young_male'; offset = 169; % cut p-wave
% patient = 'normal_male'; offset = 1;

dep = loadmat(append(DATA_PATH, 'Dataset\ECGSIM_', patient, '/ventricular_beats/beat1/user.dep'));
rep = loadmat(append(DATA_PATH, 'Dataset\ECGSIM_', patient, '/ventricular_beats/beat1/user.rep'));

A       = loadmat(append(DATA_PATH, 'Dataset\ECGSIM_', patient, '/leads/ventricles2standard_12.mat'));
% sig_ref = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/standard_12.refECG'));
% sig_ref = sig_ref(:, offset:offset+STOPTIME-1);

% [ventri_ver, ventri_tri] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.tri'));
[ventri_epi_ver, ventri_epi_tri, ventri_epi_ver_idx] = loadtri_ecgsim(append(DATA_PATH, 'Dataset\ECGSIM_', patient, '/model/ventricle_epi.tri'));
[ventri_endo_ver, ventri_endo_tri, ventri_endo_ver_idx] = loadtri_ecgsim(append(DATA_PATH, 'Dataset\ECGSIM_', patient, '/model/ventricle_endo.tri'));


CT = ones(length(ventri_epi_ver_idx) + length(ventri_endo_ver_idx),1);
for i=1:length(ventri_endo_ver_idx)
    CT(ventri_endo_ver_idx(i)) = 3;
end

N_nodes = length(dep);
TMP_matrix = zeros(N_nodes, STOPTIME/HT);
TMP_ref = zeros(N_nodes, STOPTIME/HT);
T_mapped = zeros(N_nodes, STOPTIME/HT);
ECG = zeros(9, 3, 12, STOPTIME/HT);

% treshold to calculate time points of depolarization and repolarization
threshold = -40; % [mV]

name = ["epicardium", "global", "endocardium"];

mask_ep = (CT == 1);
mask_ed = (CT == 3);

% phases currents modification (x rows for x phases)
phases_mod = [
    0 1 2;
    0.8 1 1.2;
    0.8 1 1.2
    ];

[X, Y] = ndgrid(1:size(phases_mod,1), 1:size(phases_mod,2));
combinations = [X(:), Y(:)];


RMSE = zeros(3, size(phases_mod,1), size(phases_mod,2));
CORR = zeros(3, size(phases_mod,1), size(phases_mod,2), 12);
RD = zeros(3, size(phases_mod,1), size(phases_mod,2), 12);
% matrixes to store results:
% 1. dim = epi, global, endo
% 2. dim = phase number
% 3. dim = factor number

find_fiducials = @(t, V) deal(...
    t(find(V >= threshold, 1, 'first')), ... % depolarization time
    t(find(V(find(V >= threshold, 1, 'first'):end) <= threshold, 1, 'first') + find(V >= threshold, 1, 'first') - 1) ... % repolarization time
);


% getting reference signals (without modification)
[t_tmpl_epi, V_tmpl_epi] = wrapper_TenTusscher2mod(HT, STOPTIME, 1, [1 1 1]);
[t_dep_epi, t_rep_epi] = find_fiducials(t_tmpl_epi, V_tmpl_epi);
[t_tmpl_endo, V_tmpl_endo] = wrapper_TenTusscher2mod(HT, STOPTIME, 3, [1 1 1]);
[t_dep_endo, t_rep_endo] = find_fiducials(t_tmpl_endo, V_tmpl_endo);


for i = 1:N_nodes
    if CT(i) == 1
        t_tmpl = t_tmpl_epi; V_tmpl = V_tmpl_epi;
        t_dep = t_dep_epi; t_rep = t_rep_epi;
    elseif CT(i) == 3
        t_tmpl = t_tmpl_endo; V_tmpl = V_tmpl_endo;
        t_dep = t_dep_endo; t_rep = t_rep_endo;
    else
        error(sprintf("unknown cell type: %d", CT(i)));
    end
    dep_i = dep(i);
    rep_i = rep(i);

    % reference AP are interpolated that dep and rep times fit measured
    scale = (t_rep - t_dep) / (rep_i - dep_i);
    shift = t_dep - scale * dep_i;
    t_mapped = scale * t_tmpl + shift;
    
    TMP_ref(i, :) = interp1(t_tmpl, V_tmpl, t_mapped, 'linear', 'extrap');
    T_mapped(i, :) = t_mapped(:);
end
sig_ref = A * TMP_ref;

    
for t = 1:3 % type to modification epi, global, endo
    for c=1:size(combinations,1) % combinatins of modified parameters
        phase_factor_n = combinations(c,1);
        phase_n = combinations(c,2);
    
        params = ones(3);
        params(phase_n) = phases_mod(phase_n, phase_factor_n);

        [t_tmpl_epi_mod, V_tmpl_epi_mod] = wrapper_TenTusscher2mod(HT, STOPTIME, 1, params);
        [t_dep_epi_mod, t_rep_epi_mod] = find_fiducials(t_tmpl_epi_mod, V_tmpl_epi_mod);
        [t_tmpl_endo_mod, V_tmpl_endo_mod] = wrapper_TenTusscher2mod(HT, STOPTIME, 3, params);
        [t_dep_endo_mod, t_rep_endo_mod] = find_fiducials(t_tmpl_endo_mod, V_tmpl_endo_mod);
        
        % progressbar
        printProgress((t-1) * size(combinations,1) + c, 3*size(combinations,1), sprintf("mod %s params: %.1f %.1f %.1f", name(t), params(1), params(2), params(3)))
    
        for i = 1:N_nodes
            % check if curren node's tissue is modified:
            % t==2  -> global params modifiaction
            % t==CT -> current CellType params modification
            is_modified = (t == 2) || (t == 1 && CT(i) == 1) || (t == 3 && CT(i) == 3);
            if is_modified
                % Selecting the appropriate AP template for MODIFIED based on CT
                if CT(i) == 1
                    t_tmpl = t_tmpl_epi_mod; V_tmpl = V_tmpl_epi_mod;
                elseif CT(i) == 3
                    t_tmpl = t_tmpl_endo_mod; V_tmpl = V_tmpl_endo_mod;
                else
                    error(sprintf("unknown cell type in node %d: %d", i, CT(i)));
                end
    
                t_mapped = T_mapped(i, :);        
                result_V_tmpl = interp1(t_tmpl, V_tmpl, t_mapped, 'linear', 'extrap');
            else
                % Selecting the appropriate AP template for UNMODIFIED based on CT
                result_V_tmpl = TMP_ref(i, :);
            end
            TMP_matrix(i, :) = result_V_tmpl;
        end

    
        sig_sim = A * TMP_matrix;
        for j = 1:12
            r = corrcoef(sig_ref(j,:), sig_sim(j,:));
            CORR(t, phase_factor_n, phase_n, j) = r(1,2);
            RD(t, phase_factor_n, phase_n, j) = RelativeDistance(sig_ref(j,:), sig_sim(j,:));
        end

        % store generated ECG signals
        ECG((t-1)*3 + phase_n, phase_factor_n,:,:) = sig_sim;

    end
end

%% ABSTRACT - summary of the effects of parameter modification for the selected lead
disp('Running standalone AP modification simulations for the ABSTRACT figure...');

% 1. Setup & Load Data (Independent block)
num_subjects = 1;
subject = sprintf('IKEM_Pat%03d', num_subjects);
path = 'C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data\Dataset\';

% Load geometry and signals
DATA = readGeomPeacsModelDataset(path, subject);
CineEcgResultPath = fullfile(path, subject, 'signals', sprintf('IKEM_Pat%03d.iecg', num_subjects));
EcgDataPath = fullfile(path, subject, 'signals', 'ECG_DATA');
N = 4; % RV pacing
cine_ecg = readCineEcgResults(CineEcgResultPath, 'ecgdir', EcgDataPath, 'domedian');

% Build minimal GEOM structure for initialization
GEOM = struct();
GEOM.VER = DATA.VENTR.geom.VER;
GEOM.AMA = DATA.VENTR.AMA.x99Prague;
GEOM.BSM = cine_ecg.MEDIANDATA.VENTRICULAR{N}.beats.ECGbeat;

% 2. Initialize Depolarization and Repolarization Times
disp('Initializing depolarization times (simplyfocusscan)...');
init_scan_dep = simplyfocusscan(GEOM, 1.0);
init_dep = init_scan_dep.dep;

disp('Initializing repolarization times (lineargradientscan)...');
init_scan_rep = lineargradientscan(GEOM, [-1, 0, 0], 0.5);
init_rep = init_scan_rep.rep;

% 3. Simulation Parameters
HT = 0.5;               % Time step [ms]
STOPTIME = 600;         % Simulation duration [ms]
N_nodes = length(init_dep);
threshold = -40;        % [mV] threshold for AP boundary detection
lead_to_plot = 6;       % Lead index (e.g., 6 for aVF)

% Map cell types (1 = Epicardium, 3 = Endocardium) based on endoVER array
CT = ones(N_nodes, 1);
CT(DATA.GEOM.ventr.endoVER > 0) = 3;

% Helper function to find AP boundaries (depolarization and repolarization times)
find_fiducials = @(t, V) deal(...
    t(find(V >= threshold, 1, 'first')), ...
    t(find(V(find(V >= threshold, 1, 'first'):end) <= threshold, 1, 'first') + find(V >= threshold, 1, 'first') - 1) ...
);

% 4. Generate Reference AP Templates and Mapped TMPs
[t_tmpl_epi, V_tmpl_epi] = wrapper_TenTusscher2mod(HT, STOPTIME, 1, [1 1 1]);
[t_dep_epi, t_rep_epi] = find_fiducials(t_tmpl_epi, V_tmpl_epi);

[t_tmpl_endo, V_tmpl_endo] = wrapper_TenTusscher2mod(HT, STOPTIME, 3, [1 1 1]);
[t_dep_endo, t_rep_endo] = find_fiducials(t_tmpl_endo, V_tmpl_endo);

TMP_ref = zeros(N_nodes, STOPTIME/HT);
T_mapped = zeros(N_nodes, STOPTIME/HT);

% ZDEFINIOWANA WSPÓLNA OŚ CZASU (1 x 1200)
t_standard = linspace(HT, STOPTIME, STOPTIME/HT);

for i = 1:N_nodes
    if CT(i) == 1
        t_tmpl = t_tmpl_epi; V_tmpl = V_tmpl_epi;
        t_dep = t_dep_epi; t_rep = t_rep_epi;
    else
        t_tmpl = t_tmpl_endo; V_tmpl = V_tmpl_endo;
        t_dep = t_dep_endo; t_rep = t_rep_endo;
    end

    % Correct mathematical mapping: scale the standard time axis into the template's time frame
    scale = (t_rep - t_dep) / (init_rep(i) - init_dep(i));
    shift = t_dep - scale * init_dep(i);
    t_mapped = scale * t_standard + shift;

    % Interpolate template values strictly at the standard time points
    TMP_ref(i, :) = interp1(t_tmpl, V_tmpl, t_mapped, 'linear', 'extrap');
    T_mapped(i, :) = t_mapped;
end

% Reference ECG calculation
ecg_ref_all = GEOM.AMA * TMP_ref;

% 5. Define modification variants and set up plot
phases_mod = [
    0.5 1.5; % Phase 1: negative (-) / positive (+) multiplier
    0.8 1.2; % Phase 2: negative (-) / positive (+) multiplier
    0.8 1.2  % Phase 3: negative (-) / positive (+) multiplier
];

y_labels = {'Epicardium', 'Global', 'Endocardium'};
x_labels = {'Phase 1', 'Phase 2', 'Phase 3'};

figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

% ZMIANA: Użycie cell array {} rozbija tytuł na dwie linie bez błędów parsera TeX
sgtitle({'Comparison of the influence of modification of currents of individual phases', ...
         'of action potential (AP) repolarization on the simulation of RMS of ECG signals'}, ...
        'FontSize', 16, 'FontWeight', 'bold');

time_axis = linspace(HT, STOPTIME, STOPTIME/HT);

% 6. Main simulation loop
for layer_idx = 1:3 % 1: Epicardium, 2: Global, 3: Endocardium
    for phase_idx = 1:3 % 1: Phase 1, 2: Phase 2, 3: Phase 3
        
        num = (layer_idx - 1) * 3 + phase_idx;
        subplot(3, 3, num);
        hold on; grid on;
        
        % Negative modification templates
        params_neg = [1 1 1]; params_neg(phase_idx) = phases_mod(phase_idx, 1);
        [t_n_epi, V_n_epi] = wrapper_TenTusscher2mod(HT, STOPTIME, 1, params_neg);
        [t_n_endo, V_n_endo] = wrapper_TenTusscher2mod(HT, STOPTIME, 3, params_neg);
        
        % Positive modification templates
        params_pos = [1 1 1]; params_pos(phase_idx) = phases_mod(phase_idx, 2);
        [t_p_epi, V_p_epi] = wrapper_TenTusscher2mod(HT, STOPTIME, 1, params_pos);
        [t_p_endo, V_p_endo] = wrapper_TenTusscher2mod(HT, STOPTIME, 3, params_pos);
        
        TMP_neg = TMP_ref;
        TMP_pos = TMP_ref;
        
        % Substitute modified TMPs for relevant nodes
        for i = 1:N_nodes
            is_modified = (layer_idx == 2) || (layer_idx == 1 && CT(i) == 1) || (layer_idx == 3 && CT(i) == 3);
            
            if is_modified
                t_mapped = T_mapped(i, :);
                if CT(i) == 1 % Epi
                    TMP_neg(i, :) = interp1(t_n_epi, V_n_epi, t_mapped, 'linear', 'extrap');
                    TMP_pos(i, :) = interp1(t_p_epi, V_p_epi, t_mapped, 'linear', 'extrap');
                else % Endo
                    TMP_neg(i, :) = interp1(t_n_endo, V_n_endo, t_mapped, 'linear', 'extrap');
                    TMP_pos(i, :) = interp1(t_p_endo, V_p_endo, t_mapped, 'linear', 'extrap');
                end
            end
        end
        
        % Calculate modified ECGs
        ecg_neg_all = GEOM.AMA * TMP_neg;
        ecg_pos_all = GEOM.AMA * TMP_pos;
        
        % RMS calculation
        ecg_0 = lowpassma(rms(ecg_ref_all, 1), 50);
        ecg_n = lowpassma(rms(ecg_neg_all, 1), 50);
        ecg_p = lowpassma(rms(ecg_pos_all, 1), 50);
        
        % Plotting
        h_n = plot(time_axis, ecg_n, 'b', 'LineWidth', 1.5);
        h_0 = plot(time_axis, ecg_0, 'k', 'LineWidth', 1.5);
        h_p = plot(time_axis, ecg_p, 'r', 'LineWidth', 1.5);
        
        % ZMIANA: Nowe limity osi zgodnie z prośbą
        ylim([-20, 120]);
        xlim([0, STOPTIME]);
        
        if phase_idx == 1
            ylabel('Amplitude [mV]', 'FontSize', 12, 'FontWeight', 'bold');
        else
            yticklabels({});
        end
        
        if phase_idx == 3
            text(1.08, 0.5, y_labels{layer_idx}, 'Units', 'normalized', ...
                 'Rotation', -90, 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', 'FontSize', 14);
        end
        
        if layer_idx == 1
            title(x_labels{phase_idx}, 'FontWeight', 'bold', 'FontSize', 14);
        end
        
        if layer_idx == 3
            xlabel('Time [ms]', 'FontSize', 12, 'FontWeight', 'bold');
        else
            xticklabels({});
        end
        
        if layer_idx == 3 && phase_idx == 3
            legend([h_p, h_0, h_n], '+', 'ref', '-', 'Location', 'northeast');
        end
        
    end
end

%% summary of metrics (correlation and relative distance)
% of the effects of parameter modification for the selected lead
function plot_4d_corr(corr_data, varargin)
    % plot_4d_corr Visualizes a 3x3x3x12 correlation matrix or its average.
    
    % 1. Validate Input Size
    if ~isequal(size(corr_data), [3, 3, 3, 12])
        error('Input data must be exactly 3x3x3x12.');
    end

    % 2. Parse Optional Descriptions
    p = inputParser;
    addParameter(p, 'Dim1Labels', {'Y1', 'Y2', 'Y3'}, @iscell);
    addParameter(p, 'Dim2Labels', {'X1', 'X2', 'X3'}, @iscell);
    addParameter(p, 'Dim3Labels', {'Block 1', 'Block 2', 'Block 3'}, @iscell);
    addParameter(p, 'AverageLeads', false, @islogical); % Nowa opcja
    
    default_leads = {'I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};
    addParameter(p, 'Dim4Labels', default_leads, @iscell);
    addParameter(p, 'Title', '4D Correlation Analysis', @ischar);
    
    parse(p, varargin{:});
    opts = p.Results;

    % 3. Handle Averaging logic
    if opts.AverageLeads
        % Średnia po 4. wymiarze (odprowadzenia) -> wynik 3x3x3x1
        plot_data = mean(corr_data, 4);
        num_plots = 1;
        fig_pos = [200, 200, 800, 500]; % Mniejsze okno dla pojedynczego wykresu
    else
        plot_data = corr_data;
        num_plots = 12;
        fig_pos = [100, 100, 1400, 800];
    end

    % 4. Set up the figure
    figure('Color', 'w', 'Position', fig_pos);
    
    % 5. Loop through plots
    for i = 1:num_plots
        if num_plots > 1
            subplot(3, 4, i);
        end
        hold on;

        lead_data = squeeze(plot_data(:, :, :, i));
        % Konkatenacja bloków Dim3 poziomo
        concat_data = [lead_data(:,:,1), lead_data(:,:,2), lead_data(:,:,3)];
        
        imagesc(concat_data);
        colormap(gca, parula); 
        clim([0, 1]); % Ustawione zgodnie z Twoim oryginalnym kodem (0 do 1)

        % Labels and Formatting
        yticks(1:3);
        yticklabels(opts.Dim1Labels);
        xticks([2, 5, 8]); 
        xticklabels(opts.Dim3Labels); 
        
        % Separatory między blokami
        plot([3.5, 3.5], [0.5, 3.5], 'k-', 'LineWidth', 2);
        plot([6.5, 6.5], [0.5, 3.5], 'k-', 'LineWidth', 2);
        
        if opts.AverageLeads
            title('Average of All Leads', 'FontWeight', 'bold');
        else
            title(opts.Dim4Labels{i}, 'FontWeight', 'bold');
        end
        
        axis ij; 
        axis tight;
        
        % Overlay numerical values
        for r = 1:3
            for c = 1:9
                val = concat_data(r, c);
                txt_color = 'k';
                if val > 0.6; txt_color = 'w'; end
                text(c, r, sprintf('%.2f', val), 'HorizontalAlignment', 'center', ...
                     'Color', txt_color, 'FontSize', 9, 'FontWeight', 'bold');
            end
        end
    end

    % Global Colorbar
    if opts.AverageLeads
        colorbar; % Standardowy colorbar obok wykresu
    else
        cb = colorbar('Position', [0.93 0.1 0.02 0.8]);
        cb.Label.FontSize = 12;
    end
    
    % sgtitle(opts.Title, 'FontSize', 16, 'FontWeight', 'bold');
    annotation('textbox', [0, 0.97, 1, 0.03], ...
        'String', opts.Title, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center');
end

% close all

y_labels = {'Epi', 'Global', 'Endo'};           % Dim 1
x_labels = {'-', '1', '+'};   % Dim 2
cond_labels = {'Phase 1', 'Phase 2', 'Phase 3'};                  % Dim 3

plot_4d_corr(CORR, ...
    'Dim1Labels', y_labels, ...
    'Dim2Labels', x_labels, ...
    'Dim3Labels', cond_labels, ...
    'Title', 'Correlation', ...
    'AverageLeads', true ...
    );

plot_4d_corr(RD, ...
    'Dim1Labels', y_labels, ...
    'Dim2Labels', x_labels, ...
    'Dim3Labels', cond_labels, ...
    'Title', 'RelativeDistance', ...
    'AverageLeads', true ...
    );

%% AP preview for three randomly selected nodes (interpolation validation)
% close all

figure;
hold on;
rand_nodes = randi(N_nodes, 3, 1);
t = linspace(HT,STOPTIME,STOPTIME/HT);
colors = lines(3);
for k = 1:3
    node = rand_nodes(k);
    plot(t, TMP_ref(node, :), 'Color', colors(k,:), ...
        'DisplayName', sprintf('Node %d: %s (Dep: %.1f, Rep: %.1f)', node, name(CT(node)), dep(node), rep(node)));
    plot(dep(node), -40, 'o', 'Color', colors(k,:), 'HandleVisibility', 'off');
    plot(rep(node), -40, 'x', 'Color', colors(k,:), 'HandleVisibility', 'off');
end
yline(-40, 'k--', 'threshold -40mV', 'HandleVisibility', 'off');
xlabel('t (ms)');
ylabel('V_{membrane} (mV)');
title('TMP');
legend('Location', 'best');
grid on;
hold off;


%% plot AP modifications - Aalborg 2026 presentation
% close all

y_labels = {'Epicardium', 'Mid-myocardial', 'Endocardium'};
x_labels = {'Phase 1', 'Phase 2', 'Phase 3'};

x_len = length(x_labels);
y_len = length(y_labels);

p = [1, 1, 1];
HT = 0.02;
figure

t = linspace(HT,STOPTIME,STOPTIME/HT);
params_mod = [
    0.5 1.5,
    0.8 1.2,
    0.8 1.2,
    ];

for i = 1:x_len
    for j = 1:y_len
        num = (j-1)*x_len+i;
        subplot(y_len, x_len, num)
        hold on;
        grid on;

        % CT = 1;
        % if j == 2
        %     CT = 3;
        % end
        CT = j;

        % ylim([-300, 500]);
        params = p;
        [t_tmpl_mod_0, V_tmpl_mod_0] = wrapper_TenTusscher2mod(HT, STOPTIME, j, params);
        h0 = plot( ...
            t_tmpl_mod_0, V_tmpl_mod_0, 'k', ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf("params= %f %f %f", params(1), params(2), params(3)));

        params(i) = params_mod(i, 1);
        [t_tmpl_mod_p, V_tmpl_mod_p] = wrapper_TenTusscher2mod(HT, STOPTIME, j, params);
        h1 = plot( ...
            t_tmpl_mod_p, V_tmpl_mod_p, 'r', ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf("params= %f %f %f", params(1), params(2), params(3)));

        params(i) = params_mod(i, 2);
        [t_tmpl_mod_n, V_tmpl_mod_n] = wrapper_TenTusscher2mod(HT, STOPTIME, j, params);
        h3 = plot( ...
            t_tmpl_mod_n, V_tmpl_mod_n, 'b', ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf("params= %f %f %f", params(1), params(2), params(3)));

        ylabel('Amplitude [mV]', 'FontSize', 10);
        % yticklabels({});
    
        if i == x_len
            text(1.08, 0.5, y_labels{j}, 'Units', 'normalized', ...
                     'Rotation', -90, 'HorizontalAlignment', 'center', ...
                     'FontWeight', 'bold', 'FontSize', 12);
        end
        if j == 1
            title(x_labels{i}, 'FontWeight', 'bold', 'FontSize', 12);
        end
        xlabel('Time [ms]', 'FontSize', 10);
        % xticklabels({});
        xlim([-50, 750]);

        lgd = legend([h1, h0, h3], sprintf('x %.1f', params_mod(i, 1)), 'unchanged', sprintf('x %.1f', params_mod(i, 2)));
    end
end

%% ------------------------------------------------------------------------
% -------------------- CINC2026 subjects comparison -----------------------
% -------------------------------------------------------------------------

subjects = [1,5,6,7,8];

heart_ver = {};
heart_tri = {};
heart_val = {};
thorax_ver = {};
thorax_tri = {};
thorax_val = {};
thorax_leads = {};
rd = [];

for i = 1:size(subjects,2)
    subject_num = subjects(i);
    subject = sprintf('IKEM_Pat%03d', subject_num);
    fprintf(newline+"_________ processing %s ________"+newline, subject);

    path = 'C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data\Dataset\';
    DATA = readGeomPeacsModelDataset(path, subject );    
    CineEcgResultPath = fullfile(path, subject, 'signals', sprintf('IKEM_Pat%03d.iecg', subject_num));
    EcgDataPath = fullfile(path, subject, 'signals', 'ECG_DATA');
    
    N = 4; % 3 for LV pacing or 4 for RV pacing
    cine_ecg = readCineEcgResults(CineEcgResultPath,'ecgdir', EcgDataPath, 'domedian'); 
    BSM = cine_ecg.MEDIANDATA.VENTRICULAR{N}.beats.ECGbeat;
    [TST, init_values, thorax_pot] = get_results_for_single_patient(DATA, BSM);
    thorax_ref_v = thorax_pot{1};
    thorax_sim_v = thorax_pot{2};
    diff = thorax_ref_v - thorax_sim_v;

    heart_ver = [heart_ver, DATA.VENTR.geom.VER];
    heart_tri = [heart_tri, DATA.VENTR.geom.ITRI];
    heart_val = [heart_val, init_values(2)];
    thorax_ver = [thorax_ver, DATA.GEOM.thorax.VER];
    thorax_tri = [thorax_tri, DATA.GEOM.thorax.ITRI];
    thorax_val = [thorax_val, diff];
    thorax_leads = [thorax_leads, DATA.VENTR.LEADPOS.x99Prague];
    rd = [rd, TST.rd];
end


%%
q = initQtripy();
q.reset();
q.disable_debounce();
q.background_color("white");
q.set_panels_number(1,size(subjects,2));

if true % thorax
    q.text(sprintf("ref - sim"), [0.78, 0.95]);
    for i = 1:size(subjects, 2)
        q.text(sprintf("RD=%.3f", rd(i)), [0.05, i/size(subjects, 2) - 0.1]);
        q.set_active_panel(1,i);

        offset = mean(thorax_ver{i}, 1);
        thorax_ver{i} = thorax_ver{i} - offset;
        thorax_leads{i} = thorax_leads{i} - offset;

        q.surface(thorax_ver{i}, thorax_tri{i});
        q.transparency(0.0);
        q.values(thorax_val{i});
        q.gradient_bins(10);
        q.markers(thorax_leads{i}, 'black', 10);
    end
    q.color_range(-1, 1);
end

if false % heart
    q.text(sprintf("rep time [ms]"), [0.65, 0.95]);
    for i = 1:size(subjects, 2)
        q.set_active_panel(1,i);

        offset = mean(heart_ver{i}, 1);
        heart_ver{i} = heart_ver{i} - offset;

        q.surface(heart_ver{i}, heart_tri{i});
        q.transparency(0.0);
        q.values(heart_val{i});
        q.gradient_bins(5);
        q.cmd("angle 0 110 -50");
    end

    max_in_cell = cellfun(@(x) max(x(:)), heart_val);
    max_val = max(max_in_cell);
    min_in_cell = cellfun(@(x) min(x(:)), heart_val);
    min_val = min(min_in_cell);
    q.color_range(min_val, max_val);
end



%% HELPERS

function plot_ecg_signals(PHIREF, PHIA, LAY, rd)
    figure(11); clf;
    set(gcf, 'Name', 'Measured vs Simulated ecg');
    try
        sigplot(PHIREF, 'Blue = Measured, Red = Simulated', LAY, 1.5, 'b', 1);
        hold on;
        sigplot(PHIA, '', LAY, 1.5, 'r', 1);
        title(sprintf('Blue = Measured, Red = Simulated | RD = %4.4f', rd));
    catch
        plot(PHIREF(:, 1), 'b', 'LineWidth', 1.5); hold on;
        plot(PHIA(:, 1), 'r', 'LineWidth', 1.5);
        title('Blue = Measured, Red = Simulated');
    end
    drawnow;
end

function [TST, init_values, thorax_pot] = get_results_for_single_patient(DATA, BSM)
    STOPTIME = 500;

    GEOM = struct();
    GEOM.subject = DATA.subject;
    GEOM.type = 'ventricles';
    
    lead_system = 'x99Prague'; %'x65Nijmegen' 'x12plus3leads';
    
    GEOM.VER        = DATA.VENTR.geom.VER;
    GEOM.ITRI       = DATA.VENTR.geom.ITRI;
    GEOM.AMA        = DATA.VENTR.AMA.(lead_system);
    GEOM.BSM        = BSM;
    GEOM.DIST       = DATA.VENTR.DIST3D;
    GEOM.DISTsurf   = DATA.VENTR.DISTsurf;
    GEOM.DIST2W     = DATA.VENTR.DISTanis;
    GEOM.ADJ2W      = DATA.VENTR.ADJanis;   % ventricles.adjanis; 
    GEOM.neigh      = DATA.VENTR.ADJneigh;  % ventricle.adjneigh
    GEOM.ADJ        = DATA.VENTR.ADJ3D;     % ventricle.adj3d
    GEOM.LAY        = loadmat('C:\Users\Admin\Documents\Projects\ecg_project\Scripts\Matlab\ecg_analysis\inverseArno\BEM\inverse\mla\prague99.mla');
    GEOM.RegionIdx  = DATA.GEOM.ventr.segments;
    
    num_nodes           = length(GEOM.VER);
    GEOM.TVER           = zeros(size(GEOM.BSM, 1), 3);
    GEOM.purkinjever    = zeros(num_nodes, 1);
    GEOM.Lpurkinjever   = zeros(num_nodes, 1);
    GEOM.Rpurkinjever   = zeros(num_nodes, 1);
    GEOM.Rfreewallver   = zeros(num_nodes, 1);
    GEOM.endoVER        = DATA.GEOM.ventr.endoVER;
    
    
    GEOM.pS = [];
    GEOM.anisotropyRatio = 0.5; % value only to save in logs
    
    % GEOM.specs(2) - początek analizowanego sygnału (np. początek QRS)
    % GEOM.specs(3) - koniec QRS (na tej podstawie liczone jest qrsduration)
    % GEOM.specs(4) - używane do estymacji repolaryzacji
    % GEOM.specs(5) - koniec analizowanego sygnału (koniec załamka T)
    GEOM.specs = [0, 1, 100, 200, STOPTIME];
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
    GEOM.LUT_EPI = getApLut_ResizedApd(150, 475, 1, 1);
    GEOM.LUT_ENDO = getApLut_ResizedApd(150, 475, 1, 3);
    GEOM.LUT = GEOM.LUT_EPI;

    % fociscan
    velocity = 1.0; % m/s
    init = simplyfocusscan(GEOM, velocity);
    init_dep = init.dep;
    
    % ___ LINEAR GRADIENT ___
    % rep_results = get_init_rep_lineargradient(GEOM, [-1, 0, 0], 0.5);
    % init_rep = rep_results.init_rep;

    % ___ EPI ENDO ___
    % rep_results = get_init_rep_epiendo(GEOM, init_dep);
    % init_rep = rep_results.init_rep;

    % ___ COMBINED ___
    rep_results = get_init_rep_combined(GEOM, init_dep, [-1, 0, 0], 0.5);
    init_rep = rep_results.init_rep;

    init_values = {init_dep, init_rep};
    
    INV = struct();
    INV.AMA   = GEOM.AMA;
    INV.SPECS = GEOM.SPECS;
    INV.RegionIdx = DATA.GEOM.ventr.segments;
    INV.lpass = 1; % Default lowpass filter setting
    INV.useWeighedRd = 0;
    INV.LUT = GEOM.LUT;

    % INV.LUT = fixed_rep_results.assigned_LUT; % get_init_rep_epiendo
    
    INV.BSM = GEOM.BSM;
    INV.PHIREF = INV.BSM;
    INV.normphi = norm(INV.PHIREF, 'fro');
    
    usetimes = size(INV.PHIREF, 2);
    INV.T = ones(size(INV.AMA, 2), 1) * (0 : usetimes - 1);
    
    [INV.REGOP, INV.REGOPREP] = calcREGOP(GEOM, 1);
    INV.MuValues = [0, 0]; % Dummy mu values
    
    % Pack initial times into a cell array and setup defaults for notch/ampl
    params_cell = {init_dep, init_rep};
    notchopt    = struct('pot', zeros(size(GEOM.VER,1),1));
    amplopt     = struct('pot', ones(size(GEOM.VER,1),1));
    
    % forward model
    TST = gettres_v_nparams(INV, params_cell, notchopt, amplopt);

    half_time = round(size(INV.PHIREF, 2) / 2);
    rms_ref = rms(INV.PHIREF, 1);
    [~, t_peak_offset] = max(rms_ref(half_time:end));
    t_peak = half_time + t_peak_offset - 1;
    
    % 2. Pobranie wektorów potencjałów dla znalezionego czasu
    pot_ref_99 = INV.PHIREF(:, t_peak); % Potencjały referencyjne (99 odprowadzeń)
    pot_sim = DATA.VENTR.THORAX * TST.S(:, t_peak); % Symulacja na pełnej klatce
    
    % Wykonujemy błyskawiczną estymację potencjałów serca z 99 elektrod
    % przy użyciu macierzy INV.AMA i regularyzacji Tichonowa
    lambda_vis = 1e-2; 
    pot_heart_est = (INV.AMA' * INV.AMA + lambda_vis * eye(size(INV.AMA, 2))) \ (INV.AMA' * pot_ref_99);
    
    % Rzutujemy uzyskane potencjały serca na PEŁNĄ siatkę klatki piersiowej
    pot_ref_full = DATA.VENTR.THORAX * pot_heart_est;
    % ---------------------------------------
    
    thorax_pot = {pot_sim, pot_ref_full};
    
end

% function forward_ecg_with_modified_params(DATA, params)
%     CT = 1;
%     HT = 0.01;
%     STOPTIME = 800;
%     ECG         = DATA.VENTR.SIGNALS.ECG.s1(:,800:800+STOPTIME-1);
% 
%     [t, V] = wrapper_TenTusscher2mod(HT, STOPTIME, CT, params);
% 
% end