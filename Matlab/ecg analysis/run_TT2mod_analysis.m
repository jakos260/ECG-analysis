close all

% time setup
HT = 0.1;                       % timestep (ms)
STOPTIME = 500; % Duration of the simulation (ms)

% data read
patient = 'normal_young_male'; offset = 169; % cut p-wave
% patient = 'normal_male'; offset = 1;

dep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.dep'));
rep = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ventricular_beats/beat1/user.rep'));

A       = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricles2standard_12.mat'));
% sig_ref = loadmat(append(DATA_PATH, 'ECGsim_data/', patient, '/ecgs/standard_12.refECG'));
% sig_ref = sig_ref(:, offset:offset+STOPTIME-1);

% [ventri_ver, ventri_tri] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.tri'));
[ventri_epi_ver, ventri_epi_tri, ventri_epi_ver_idx] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle_epi.tri'));
[ventri_endo_ver, ventri_endo_tri, ventri_endo_ver_idx] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle_endo.tri'));


CT = ones(length(ventri_epi_ver_idx) + length(ventri_endo_ver_idx),1);
for i=1:length(ventri_endo_ver_idx)
    CT(ventri_endo_ver_idx(i)) = 3;
end

N_nodes = length(dep);
TMP_matrix = zeros(N_nodes, STOPTIME/HT);
TMP_ref = zeros(N_nodes, STOPTIME/HT);
T_mapped = zeros(N_nodes, STOPTIME/HT);
ECG = zeros(9, 3, 12, STOPTIME/HT);

threshold = -40; % [mV]

name = ["epicardium", "global", "endocardium"];

mask_ep = (CT == 1);
mask_ed = (CT == 3);

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
% 1.d = epi, global, endo
% 2.d = phase number
% 3.d = factor number

find_fiducials = @(t, V) deal(...
    t(find(V >= threshold, 1, 'first')), ... % depolarization time
    t(find(V(find(V >= threshold, 1, 'first'):end) <= threshold, 1, 'first') + find(V >= threshold, 1, 'first') - 1) ... % repolarization time
);



[t_tmpl_epi, V_tmpl_epi] = TenTusscher2mod(HT, STOPTIME, 1, [1 1 1]);
[t_dep_epi, t_rep_epi] = find_fiducials(t_tmpl_epi, V_tmpl_epi);
[t_tmpl_endo, V_tmpl_endo] = TenTusscher2mod(HT, STOPTIME, 3, [1 1 1]);
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

    scale = (t_rep - t_dep) / (rep_i - dep_i);
    shift = t_dep - scale * dep_i;
    t_mapped = scale * t_tmpl + shift;

    TMP_ref(i, :) = interp1(t_tmpl, V_tmpl, t_mapped, 'linear', 'extrap');
    T_mapped(i, :) = t_mapped(:);
end
sig_ref = A * TMP_ref;

    
for t = 1:3 % type to modification epi, global, endo
    for c=1:size(combinations,1)
        phase_factor_n = combinations(c,1);
        phase_n = combinations(c,2);
    
        params = ones(3);
        params(phase_n) = phases_mod(phase_n, phase_factor_n);

        [t_tmpl_epi_mod, V_tmpl_epi_mod] = TenTusscher2mod(HT, STOPTIME, 1, params);
        [t_dep_epi_mod, t_rep_epi_mod] = find_fiducials(t_tmpl_epi_mod, V_tmpl_epi_mod);
        [t_tmpl_endo_mod, V_tmpl_endo_mod] = TenTusscher2mod(HT, STOPTIME, 3, params);
        [t_dep_endo_mod, t_rep_endo_mod] = find_fiducials(t_tmpl_endo_mod, V_tmpl_endo_mod);
        
        printProgress((t-1) * size(combinations,1) + c, 3*size(combinations,1), sprintf("mod %s params: %.1f %.1f %.1f", name(t), params(1), params(2), params(3)))
        % fprintf("%d %d %d\n", t, phase_n, phase_factor_n);
    
        for i = 1:N_nodes
            is_modified = (t == 2) || (t == 1 && CT(i) == 1) || (t == 3 && CT(i) == 3);
            if is_modified
                % Wybór odpowiedniego template'u dla MODYFIKOWANYCH na podstawie CT
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
                % Wybór odpowiedniego template'u dla NIEMODYFIKOWANYCH na podstawie CT
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
        
        ECG((t-1)*3 + phase_n, phase_factor_n,:,:) = sig_sim;

    end
end
%%
close all

for i = 1:9
    num_to_plot = i;
    figure(num_to_plot)

    ecg_0 = lowpassma(sig_ref, 50);
    ecg_1 = lowpassma(squeeze(ECG(num_to_plot, 1, :, :)), 50);
    % ecg_2 = lowpassma(squeeze(ECG(num_to_plot, 2, :, :)), 50);
    ecg_3 = lowpassma(squeeze(ECG(num_to_plot, 3, :, :)), 50);

    leadv12( ...
        ecg_0, ecg_1, ecg_3, ...
        'paperspeed', 20, ...
        'amplification', 100 ...
        );

    n = sprintf("mod %s params of phase %d", name(combinations(i, 2)), combinations(i, 1));
    printProgress(i, 9, "plotting " + n);
    sgtitle(n)
    hold on;
    h2 = plot(NaN, NaN, 'r', 'LineWidth', 1.5);
    h1 = plot(NaN, NaN, 'k', 'LineWidth', 1.5);
    h3 = plot(NaN, NaN, 'b', 'LineWidth', 1.5);
    lgd = legend([h2, h1, h3], '+', 'ref', '-');
end

%%
close all

lead_to_plot = 6;
standard12_names = ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"];
y_labels = {'Epicardium', 'Global', 'Endocardium'};
x_labels = {'Phase 1', 'Phase 2', 'Phase 3'};

figure("Name", "Abstract figure_1")
sgtitle( ...
    sprintf("Comparison of the influence of modification of currents" + ...
    " of individual phases\nof action potential (AP) repolarization" + ...
    " on the simulation of ECG signals (lead %s)", ...
    standard12_names(lead_to_plot)), 'FontSize', 16, 'FontWeight', 'bold');

t = linspace(HT,STOPTIME,STOPTIME/HT);
for i = 1:3
    for j = 1:3
        num = (i-1)*3+j;
        subplot(3, 3, num)
        hold on;
        grid on;

        ylim([-300, 500]);

        ecg_0 = lowpassma(sig_ref(lead_to_plot, :), 50);
        ecg_1 = lowpassma(squeeze(ECG(num, 1, lead_to_plot, :)), 50);
        ecg_3 = lowpassma(squeeze(ECG(num, 3, lead_to_plot, :)), 50);
        h0 = plot(t, ecg_0, 'k', 'LineWidth', 1.5);
        h1 = plot(t, ecg_1, 'r', 'LineWidth', 1.5);
        h3 = plot(t, ecg_3, 'b', 'LineWidth', 1.5);

        if j == 1
            ylabel('Amplitude [mV]', 'FontSize', 10);
        else
            yticklabels({});
        end

        if j == 3
            text(1.08, 0.5, y_labels{i}, 'Units', 'normalized', ...
                 'Rotation', -90, 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', 'FontSize', 12);
        end
        
        if i == 1
            title(x_labels{j}, 'FontWeight', 'bold', 'FontSize', 12);
        end
        
        if i == 3
            xlabel('Time [ms]', 'FontSize', 10);
        else
            xticklabels({});
        end
    end
end
lgd = legend([h1, h0, h3], '+', 'ref', '-');


%%
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

close all

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

%%
close all

figure;
hold on;
rand_nodes = randi(N_nodes, 3, 1);
t = linspace(HT,STOPTIME,STOPTIME/HT);
colors = lines(3);
for k = 1:3
    node = rand_nodes(k);
    plot(t, TMP_ref(node, :), 'Color', colors(k,:), ...
        'DisplayName', sprintf('Node %d: %s (Dep: %.1f, Rep: %.1f)', node, name(CT(node)), dep(node), rep(node)));
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
