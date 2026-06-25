num_subjects = 1; 

subject = sprintf('IKEM_Pat%03d', num_subjects);
path = 'C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data\Dataset\';
DATA = readGeomPeacsModelDataset(path, subject);

%%
STOPTIME = 500;

GEOM = struct();
GEOM.subject = subject;
GEOM.type = 'ventricles';

lead_system = 'x99Prague'; %'x65Nijmegen' 'x12plus3leads';
ref_signal = 'x3_2019_11_15_12_36_17';

beat_no = 3;
beat_time =  DATA.VENTR.SIGNALS.INFO.(ref_signal).beats(beat_no);
beat_cut = [beat_time-60 beat_time+500];

GEOM.VER        = DATA.VENTR.geom.VER;
GEOM.ITRI       = DATA.VENTR.geom.ITRI;
GEOM.AMA        = DATA.VENTR.AMA.(lead_system);
GEOM.BSM        = DATA.VENTR.SIGNALS.BSM.(ref_signal)(:,beat_cut(1):beat_cut(2));
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

%% Get action potential LUT
load('inverseMy/TmpLut_niceApd.mat');
% LUT = getTmpLut_niceApd(200, 460, 1);
% LUT = getTmpLut_HRmod(50, 200, 1);
GEOM.LUT = LUT;

%%
close all
figure(4)
plot(GEOM.BSM(1,:));

%% multifoci skan

disp('multifociscan');
[measinit.foci, measinit.dep, measinit.outp] = multifociscan(GEOM, 1, 0);
measinit.rep = initRep(GEOM, measinit.dep);

init_dep = measinit.dep;
t_wave_peak = mean(measinit.rep);
[val, qrs_peak] = max(rms(abs(GEOM.BSM(:,1:STOPTIME))));

apd = t_wave_peak - qrs_peak;
alpha = 0.7;

init_rep = init_dep * alpha + apd;

%% inverse modeling
mudep = 1e-5;
murep = 1e-5;

meas = my_inversefunc(GEOM, init_dep, init_rep, mudep, murep);

%%

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
q.values(meas.depfinal - substract_from_dep_values);
q.gradient_bins(gradient_bins);
q.text(sprintf("dep av=%.0f[ms]", mean(meas.depfinal)), [0.1, 0.95]);

q.set_active_panel(1, 2);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(init_dep - substract_from_dep_values);
q.gradient_bins(gradient_bins);

% simulated repolarisation time
q.set_active_panel(2, 1);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(meas.repfinal - substract_from_rep_values);
q.gradient_bins(gradient_bins);
q.text(sprintf("rep av=%.0f[ms]", mean(meas.repfinal)), [0.6, 0.95]);

q.set_active_panel(2, 2);
q.surface(GEOM.VER, GEOM.ITRI);
q.transparency(transparency);
q.values(init_rep - substract_from_rep_values);
q.gradient_bins(gradient_bins);

q.color_range(0, 350);

%%
maxAmpl = round(max(max(abs(GEOM.BSM))));
figure(99);clf; sigplot(GEOM.BSM,'',GEOM.LAY,1.3/maxAmpl,'b',1,0); 
hold on; 
sigplot(meas.simulated_ecg,'',GEOM.LAY,1.3/maxAmpl,'r',1,0);	

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