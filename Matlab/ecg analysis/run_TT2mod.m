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
patient = 'normal_young_male'; offset = 169; % cut p-wave
% patient = 'normal_male'; offset = 1;

dep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.dep'));
rep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.rep'));

A       = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2BSM_(nijmegen_64).mat'));
sig_ref = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/BSM_(nijmegen_64).refECG'));
sig_ref = sig_ref(:, offset:end);

% [ventri_ver, ventri_tri] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.tri'));
[ventri_epi_ver, ventri_epi_tri, ventri_epi_ver_idx] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle_epi.tri'));
[ventri_endo_ver, ventri_endo_tri, ventri_endo_ver_idx] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle_endo.tri'));

CT = ones(length(ventri_epi_ver_idx) + length(ventri_endo_ver_idx),1);
for i=1:length(ventri_endo_ver_idx)
    CT(ventri_endo_ver_idx(i)) = 3;
end
% time setup
HT = 0.1;                       % timestep (ms)
STOPTIME = 500; % Duration of the simulation (ms)

N_nodes = length(dep);
TMP_matrix = zeros(N_nodes, STOPTIME/HT);

phases_mod = [1, 1, 1];

% --- Template Epicardial (CellType = 1) ---
[t_tmpl_epi, V_tmpl_epi] = TenTusscher2mod(HT, STOPTIME, 1, phases_mod);

% --- Template Endocardial (CellType = 3) ---
[t_tmpl_endo, V_tmpl_endo] = TenTusscher2mod(HT, STOPTIME, 3, phases_mod);

% reference points with treshold [mV]
threshold = -40;

find_fiducials = @(t, V) deal(...
    t(find(V >= threshold, 1, 'first')), ... % depolarization time
    t(find(V(find(V >= threshold, 1, 'first'):end) <= threshold, 1, 'first') + find(V >= threshold, 1, 'first') - 1) ... % repolarization time
);

[t_dep_epi, t_rep_epi] = find_fiducials(t_tmpl_epi, V_tmpl_epi);
[t_dep_endo, t_rep_endo] = find_fiducials(t_tmpl_endo, V_tmpl_endo);

i = 1;

scale = (t_rep_epi - t_dep_epi) / (rep(i) - dep(i));
shift = t_dep_epi - scale * dep(i);
time_mapped = scale * t_tmpl_epi + shift;

resampled_V_tmpl_epi = interp1(t_tmpl_epi, V_tmpl_epi, time_mapped, 'linear', 'extrap');

figure;
hold on;
plot(t_tmpl_epi, V_tmpl_epi, 'k');
plot(t_tmpl_epi, resampled_V_tmpl_epi, 'g');
yline(threshold, 'k--', sprintf("%d [mv]", threshold), 'HandleVisibility', 'off');
xline(dep(i), 'r');
xline(rep(i), 'b');
xlabel('Czas (ms)');
ylabel('V_{membrane} (mV)');
title('Zmapowane przebiegi TMP dla wybranych węzłów');
grid on;
hold off;

%%
for i = 1:N_nodes
    % Wybór odpowiedniego template'u na podstawie isEpi
    if CT(i) == 1
        t_tmpl = t_tmpl_epi; V_tmpl = V_tmpl_epi;
        t_dep = t_dep_epi; t_rep = t_rep_epi;
    elseif CT(i) == 3
        t_tmpl = t_tmpl_endo; V_tmpl = V_tmpl_endo;
        t_dep = t_dep_endo; t_rep = t_rep_endo;
    end
    
    dep_i = dep(i);
    rep_i = rep(i);
    
    scale = (t_rep - t_dep) / (rep_i - dep_i);
    shift = t_dep - scale * dep_i;
    t_mapped = scale * t_tmpl + shift;
    
    resampled_V_tmpl = interp1(t_tmpl, V_tmpl, t_mapped, 'linear', 'extrap');

    TMP_matrix(i, :) = resampled_V_tmpl;
end

figure;
hold on;
rand_nodes = randi(N_nodes, 3, 1);
colors = lines(3);
for k = 1:3
    node = rand_nodes(k);
    plot(t_tmpl_epi, TMP_matrix(node, :), 'Color', colors(k,:), ...
        'DisplayName', sprintf('Node %d (CT: %d, Dep: %.1f, Rep: %.1f)', node, CT(node), dep(node), rep(node)));
    % Zaznaczenie punktów depolaryzacji i repolaryzacji
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

