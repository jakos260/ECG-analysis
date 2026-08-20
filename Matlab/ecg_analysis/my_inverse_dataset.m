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
GEOM.subject = subject;
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
GEOM.LUT = getApLut_ResizedApd(150, 475, 1);


maxAmpl = round(max(max(abs(GEOM.BSM))));
figure(99);clf; sigplot(GEOM.BSM,'',GEOM.LAY,1.3/maxAmpl,'b',1,0); 


%% multifoci skan

% disp('multifociscan');
% % [measinit.foci, measinit.dep, measinit.outp] = multifociscan(GEOM, 1, 0);
% measinit.dep = DATA.VENTR.HEARTDIST(:,1216);
% measinit.rep = initRep(GEOM, measinit.dep);
% 
% init_dep = measinit.dep;
% % init_rep = measinit.rep;
% t_wave_peak = mean(measinit.rep);
% [val, qrs_peak] = max(rms(abs(GEOM.BSM)));
% 
% apd = t_wave_peak - qrs_peak;
% alpha = 0.8;
% 
% init_rep = init_dep * alpha + apd;

velocity = 1; % m/s
init = simplyfocusscan(GEOM, velocity);
init_dep = init.dep;
init_rep = init.rep;

%% inverse modeling
mudep = 1e-8;
murep = 2e-7;

meas = my_inversefunc(GEOM, init_dep, init_rep, mudep, murep, true);
%%
mudep = 1e-3;
murep = 2e-4;

param_configs(1).name = 'Depolarization';
param_configs(1).initial_value = init_dep;
param_configs(1).mu = mudep;
param_configs(1).bounds = [0, 150];

param_configs(2).name = 'Repolarization';
param_configs(2).initial_value = init_rep;
param_configs(2).mu = murep;
param_configs(2).bounds = [150, 400];

settings = struct();
settings.maxiter = 4;
settings.plot_en = true;
settings.lpass = 1;
settings.mode = 4;
settings.useWeighedRd = 0;
settings.MINRD = 0.18;
settings.stopcrit = 1e-4;

opt_func = @gettres_v_nparams;

inv_model   = INVERSE_MODEL(GEOM, param_configs, settings, opt_func);
results     = inv_model.run_optimization();


%% QTRIPLOT RESULTS
[dep_final, rep_final] = results.final_params{:};

transparency = 0.0;
gradient_bins = 10;
substract_from_dep_values = 0; %min(meas.depfinal);
substract_from_rep_values = 0; %min(meas.repfinal);
q = initQtripy();

q.reset();
q.disable_debounce();
q.set_panels_number(2,2);
q.background_color("white");

% simulated depolarisation time
q.set_active_panel(1, 1);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(dep_final - substract_from_dep_values);
q.gradient_bins(gradient_bins);
q.text(sprintf("dep av=%.0f[ms]", mean(dep_final)), [0.1, 0.95]);

q.set_active_panel(1, 2);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(init_dep - substract_from_dep_values);
q.gradient_bins(gradient_bins);

% simulated repolarisation time
q.set_active_panel(2, 1);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(rep_final - substract_from_rep_values);
q.gradient_bins(gradient_bins);
q.text(sprintf("rep av=%.0f[ms]", mean(rep_final)), [0.6, 0.95]);

q.set_active_panel(2, 2);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(init_rep - substract_from_rep_values);
q.gradient_bins(gradient_bins);

max_color_range = max(max(rep_final, init_rep));
q.color_range(0, max_color_range);

%%
N = 1;
offset = 1;
% Smooth the signals (using a Gaussian filter with a small window size)
window_size = 5;

% Note: If offset:offset+400 contains 401 samples, ensure simulated_ecg indices match (e.g., 1:401).
ecg_ref = smoothdata(GEOM.BSM(N, offset:offset+400), 'gaussian', window_size);
ecg_sim = smoothdata(meas.simulated_ecg(N, 1:400), 'gaussian', window_size); 
% ecg_sim_old = smoothdata(meas_old.simulated_ecg(N, 1:400), 'gaussian', window_size);

figure(67);
clf; % Clear the figure to prevent overlapping when re-running
set(gcf, 'Color', 'w');

% First subplot: New simulation (LUT based)
% subplot(2,1,1);
yline(0, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); % Add gray zero-line first so it stays in the background
hold on;
plot(ecg_ref, 'k', 'LineWidth', 2.0, 'DisplayName', 'ref'); % Thick black line for reference
plot(ecg_sim, 'r', 'LineWidth', 2.0, 'DisplayName', 'sim'); % Thick red line for simulation
legend();

% Format title to include the Relative Distance (RD) value
title(sprintf('Reference vs New LUT Optimization (RD = %.3f)', meas.rdfinal), 'FontSize', 12);
axis off; % Remove all axes, borders, and ticks

% Second subplot: Old simulation
% subplot(2,1,2);
% yline(0, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); % Add gray zero-line
% hold on;
% plot(ecg_ref, 'k', 'LineWidth', 2.0);
% plot(ecg_sim_old, 'r', 'LineWidth', 2.0);

% Format title to include the Relative Distance (RD) value from the old model
title(sprintf('Reference vs Old Optimization (RD = %.3f)', meas.rdfinal), 'FontSize', 12);
axis off; % Remove all axes, borders, and ticks