close all;

patient = 'normal_young_male'; offset = 169;
% patient = 'normal_male'; offset = 1;

% --- DATA LOADING ---
dep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.dep'));
rep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.rep'));

% A       = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2BSM_(nijmegen_64).mat'));
% sig_ref = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/BSM_(nijmegen_64).refECG'));
% sig_ref = sig_ref(:, offset:end);

A       = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2standard_12.mat'));
sig_ref = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/standard_12.refECG'));
sig_ref = sig_ref(:, offset:end);

[ventri_ver, ventri_tri] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.tri'));


% --- PARAMETERS ---
params = [0.995, 0, 0.3, 0.13, 0.1];
mode = getsMode.Two_dv_Spline;
L = 500;                            % simulation time [ms]
T = ones(length(dep),1)*(1:L);

node_idx = 20;
lead_idx = 8;
lead_names = ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"];

% --- S GENERATING ---
S = gets(T, dep, rep, params, mode);
sig_sim = A * S;
rep_sim = first_crossing_thr(S,0.5,false);


% --- DATA PLOT ---
fig = figure('Name', 'TMP demo', 'Position', [100, 100, 1200, 700]);

subplot(1, 2, 1);
plot(1:L, S(node_idx, :), 'k', 'LineWidth', 1.5);
title(sprintf('TMP - node %d', node_idx));
xlabel('time [ms]');
ylabel('tmp values');
ylim([-0.2, 1.2]);
xline(dep(node_idx), 'r--', 'Depolarization');
xline(rep(node_idx), 'b--', 'Repolarization');
grid on;

subplot(1, 2, 2);
plot(1:L, sig_sim(lead_idx, 1:L), 'g', 'LineWidth', 1.5);
hold on
plot(1:L, sig_ref(lead_idx, 1:L), 'k', 'LineWidth', 0.75);
% yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); 
title(sprintf('sim&ref ECG - %s', lead_names(lead_idx)));
xlabel('time [ms]');
ylabel('Ampl [mV]');
legend("ECG sim", "ECG ref");
xline(dep(node_idx), 'r--', 'HandleVisibility', 'off');
xline(rep(node_idx), 'b--', 'HandleVisibility', 'off');
grid on;

% --- 3D VISUALIZATION ---
q = initQtripy();
q.reset();
q.disable_debounce();
q.set_panels_number(2, 1);

% reference repolarisation time
q.set_active_panel(2, 1);
q.marker(ventri_ver(node_idx, :), 'red', 5);
q.surface(ventri_ver, ventri_tri);
q.transparency(0.3);
q.values(rep);
q.gradient_bins(10);
q.text("reference rep", [0.6, 0.85]);

% simulated repolarisation time
q.set_active_panel(1, 1);
q.marker(ventri_ver(node_idx, :), 'black', 5);
q.surface(ventri_ver, ventri_tri);
q.transparency(0.3);
q.values(rep_sim);
q.gradient_bins(10);
q.text("simulated rep", [0.1, 0.85]);