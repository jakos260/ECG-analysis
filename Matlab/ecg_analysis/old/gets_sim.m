clearvars
close all force
delete(findall(groot,'Type','uifigure'))
drawnow


q = initQtripy();
% patient = 'normal_male'; offset = 1;
patient = 'normal_young_male'; offset = 169;


dst     = loadmat(append('data/ExportData/', patient, '/model/ventricle.dst3d'));
V       = loadmat(append('data/ExportData/', patient, '/model/ventricles2Thorax.mat'));
A64     = loadmat(append('data/ExportData/', patient, '/model/ventricles2BSM_(nijmegen_64).mat'));
A12     = loadmat(append('data/ExportData/', patient, '/model/ventricles2standard_12.mat'));
BSM_ref = loadmat(append('data/ExportData/', patient, '/ecgs/BSM_(nijmegen_64).refECG'));
BSM_ref = BSM_ref(:, offset:end);
BSM_elc = loadmat(append('data/ExportData/', patient, '/ecgs/BSM_(nijmegen_64).elec'));
STD_ref = loadmat(append('data/ExportData/', patient, '/ecgs/standard_12.refECG'));
STD_ref = STD_ref(:, offset:end);
STD_elc = loadmat(append('data/ExportData/', patient, '/ecgs/standard_12.elec'));
dep     = loadmat(append('data/ExportData/', patient, '/ventricular_beats/beat1/user.dep'));
rep     = loadmat(append('data/ExportData/', patient, '/ventricular_beats/beat1/user.rep'));
L = size(BSM_ref, 2);

[ventri_ver, ventri_tri] = loadtri_ecgsim(append('data/ExportData/', patient, '/model/ventricle.tri'));
[thorax_ver, thorax_tri] = loadtri_ecgsim(append('data/ExportData/', patient, '/model/thorax.tri'));
    

% Create UI figure
fig = uifigure('Position',[100 100 600 1000]);
ax1 = uiaxes(fig,'Position',[50 170 500 200]);
ax1.NextPlot = 'add';

ax2 = uiaxes(fig,'Position',[50 420 500 250]);
hold(ax2,'on')

ax3 = uiaxes(fig,'Position',[50 720 500 250]);
hold(ax3,'on')

% --- Isolines (zero-level reference lines) ---
fig.UserData.iso2 = yline(ax2, 0, '--', ...
    'Color', [0.5 0.5 0.5], ...
    'LineWidth', 1);

fig.UserData.iso3 = yline(ax3, 0, '--', ...
    'Color', [0.5 0.5 0.5], ...
    'LineWidth', 1);

fig.UserData.time_marker = xline(ax1, 300, 'r');

[~, ~, dd_items] = get_lead_matrix();
dd = uidropdown(fig, ...
    'Position',[70 970 80 25], ...
    'Items', dd_items, ...
    'Value', 'V2', ...
    'ValueChangedFcn', @(btn,event) dropdownCallback(fig, btn, q));

% btn = uibutton(fig, 'push', ...
%     'Text', 'QTriplot', ...
%     'Position', [450 900 90 25], ...
%     'ButtonPushedFcn', @(src,~) btnCallback(fig, q));

fig.UserData.corrText2 = text(ax2, ...
    0.98, 0.02, '', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','right');

fig.UserData.corrText3 = text(ax3, ...
    0.98, 0.02, '', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','right');


fig.UserData.b = 0.5;
fig.UserData.c = 0.5;
fig.UserData.d = 0.5;
fig.UserData.number_to_plot = 20;
fig.UserData.qtrtime = 300;
fig.UserData.lead_number = 8;
fig.UserData.A64 = A64;
fig.UserData.A12 = A12;
fig.UserData.BSM_ref = BSM_ref;
fig.UserData.BSM_elc = BSM_elc;
fig.UserData.STD_ref = STD_ref;
fig.UserData.STD_elc = STD_elc;
fig.UserData.dep = dep;
fig.UserData.rep = rep;
fig.UserData.t_sim = 0;
fig.UserData.t_ref = 0;
fig.UserData.L = L;




fig.UserData.tmp     = plot(ax1,nan,nan,'k');
hold(ax2,'on')
fig.UserData.sig_ref = plot(ax2,nan,nan,'b','LineWidth',1);
fig.UserData.sig_sim = plot(ax2,nan,nan,'k','LineWidth',1.5);
hold(ax2,'off')
hold(ax2,'on')
fig.UserData.ecg_ref = plot(ax3,nan,nan,'r','LineWidth',1);
fig.UserData.ecg_sim = plot(ax3,nan,nan,'k','LineWidth',1.5);
hold(ax2,'off')


sld1 = uislider(fig,...
    'Position',[100 700 400 3],...
    'Limits',[1 size(A64, 1)],...
    'Value',fig.UserData.number_to_plot);

    % 'MajorTicks',0:20:round(mod(length(dst), 20) * 20),...
    % 'MinorTicks',[],...
sld2 = uislider(fig,...
    'Position',[100 150 400 3],...
    'Limits',[0 1],...
    'Value',fig.UserData.b);
sld3 = uislider(fig,...
    'Position',[100 100 400 3],...
    'Limits',[0 1],...
    'Value',fig.UserData.c);
sld4 = uislider(fig,...
    'Position',[100 50 400 3],...
    'Limits',[0 1],...
    'Value',fig.UserData.d);

sld5 = uislider(fig,...
    'Position',[100 400 400 3],...
    'Limits',[1 L],...
    'Value',fig.UserData.qtrtime);

addSlider(fig,sld1, q,'number_to_plot','BSM(64) idx','%d',  true);
addSlider(fig,sld2, q,'b','p0 (dep incl)','%.3f',          false);
addSlider(fig,sld3, q,'c','p1 (st  decl)','%.3f',          false);
addSlider(fig,sld4, q,'d','p2 (rep decl)','%.3f',          false);
addSlider(fig,sld5, q,'qtrtime','time to model','%d',       true);



[~, min_idx] = min(fig.UserData.dep);
focus = ventri_ver(min_idx, :);

% --- Initial draw ---
sliderCallback( ...
    fig, ...
    'number_to_plot', ...
    q, ...
    fig.UserData.number_to_plot, ...
    [], ...          % label not needed on init
    '%d', ...
    true);



% QTriplot init
q.reset();
q.disable_debounce();
q.set_panels_number(2, 1);

% reference repolarisation time
q.set_active_panel(2, 1);
q.marker(focus, 'red', 5);
q.surface(ventri_ver, ventri_tri);
q.transparency(0.3);
q.values(rep);
q.gradient_bins(10);
q.text("reference rep", [0.6, 0.85]);

% simulated repolarisation time
q.set_active_panel(1, 1);
q.marker(focus, 'black', 5);
q.surface(ventri_ver, ventri_tri);
q.transparency(0.3);
q.values(fig.UserData.t_sim);
q.gradient_bins(10);
q.text("simulated rep", [0.1, 0.85]);

q.enable_debounce(0.02)
% q.background_color('white')

function addSlider(fig, sld, q, param, name, format, isInteger)

    % Name label
    uilabel(fig,...
        'Text',name,...
        'Position',[20 sld.Position(2)-8 80 22]);

    % Value label
    lblVal = uilabel(fig,...
        'Text',sprintf(format,sld.Value),...
        'Position',[520 sld.Position(2)-8 70 22],...
        'HorizontalAlignment','right');

    % Single callback
    sld.ValueChangingFcn = @(s,e) sliderCallback(fig,param,q,e.Value,lblVal,format,isInteger);
end

function sliderCallback(fig,param,q,val,lbl,format,isInteger)

    if isInteger
        val = round(val);
    end
    if ~isempty(lbl)
        lbl.Text = sprintf(format,val);
    end

    fig.UserData.(param) = val;

    ud = fig.UserData;
    idx = ud.number_to_plot;
    
    % --- Compute signals ---
    [tmp_all, sig_sim_all] = get_h( ...
        ud.b, ud.c, ud.d, ...
        ud.dep, ud.rep, ud.L, ud.A64);

    [~, std_sim_all] = get_h( ...
        ud.b, ud.c, ud.d, ...
        ud.dep, ud.rep, ud.L, ud.A12);
    
    tmp_sim = tmp_all(idx,:);
    sig_sim = sig_sim_all(idx,:);
    sig_ref = ud.BSM_ref(idx,:);

    lead = ud.lead_number;
    ecg12_sim = std_sim_all;
    ecg12_ref = ud.STD_ref;
    sel_ecg12_sim = ecg12_sim(lead, :);
    sel_ecg12_ref = ecg12_ref(lead, :);
    
    % --- Update plots ---
    ud.tmp.XData = 1:length(tmp_sim);
    ud.tmp.YData = tmp_sim;
    
    ud.sig_sim.XData = 1:length(sig_sim);
    ud.sig_sim.YData = sig_sim;
    
    ud.sig_ref.XData = 1:length(sig_ref);
    ud.sig_ref.YData = sig_ref;

    ud.ecg_ref.XData = 1:ud.L;
    ud.ecg_ref.YData = sel_ecg12_ref;
    ud.ecg_sim.XData = 1:ud.L;
    ud.ecg_sim.YData = sel_ecg12_sim;

    ud.time_marker.Value = ud.qtrtime;
    
    % --- Correlation coefficients ---
    C2_all = corrcoef(ud.BSM_ref(:), sig_sim_all(:));
    corr2_global = C2_all(1,2);
    C2_loc = corrcoef(sig_ref, sig_sim);
    corr2_local = C2_loc(1,2);

    C3_all = corrcoef(ecg12_ref, ecg12_sim);
    corr3_global = C3_all(1,2);
    C3_loc = corrcoef(sel_ecg12_ref, sel_ecg12_sim);
    corr3_local = C3_loc(1,2);
    
    % --- Update text ---
    ud.corrText2.String = sprintf( ...
        'corr(full) = %.3f\ncorr(local) = %.3f', ...
        corr2_global, corr2_local);
    ud.corrText3.String = sprintf( ...
        'corr(full) = %.3f\ncorr(local) = %.3f', ...
        corr3_global, corr3_local);

    
    % qtriplot update
    [sim_dep, sim_rep] = getTimesFromS(tmp_all);
    ud.t_sim = sim_rep; % sim_rep, tmp_all(:,ud.qtrtime)
    q.values(ud.t_sim)

    % assign modified UserData to figure and draw
    fig.UserData = ud;
    drawnow limitrate

end

function dropdownCallback(fig, src, q)

    % Get selected lead index
    lead_idx = src.ValueIndex;

    % Store it
    fig.UserData.lead_number = lead_idx;

    % Trigger redraw using existing logic
    sliderCallback( ...
        fig, ...
        'number_to_plot', ...
        q, ...
        fig.UserData.number_to_plot, ...
        [], ...
        '%d', ...
        true);
end

function btnCallback(fig, q)
    q.values(fig.UserData.t_sim)
end

function [tmp, sig_sim] = get_h(b, c, d, dep, rep, L, A)
    maxt = L;
    T = ones(length(dep),1)*(1:maxt);
    p = [b c d]; % [-0.025, 0.2];
    
    S = twave_generator_full(T, dep, rep, p);
    tmp = S;
    sig_sim = A * S;
end

function [tmp, sig_sim] = get_h_old(b, c, d, dst, L, A64)
    v = a;          % m/s
    dst_mm = dst;     % ECGSIM stores mm
    time_ms = (dst_mm / 1000) / v * 1000;
    root_nodes = [b];
    
    % dep compute (Dijkstra)
    time_ms(time_ms == 0) = Inf;
    
    G = graph(time_ms,'upper');  % undirected graph
    
    dep = distances(G, root_nodes);
    dep = min(dep,[],1)';        % earliest arrival
    dep = dep - min(dep);        % start at 0 ms
    
    % disp([min(dep), max(dep)]);
    
    rep = 350 + mean(dep) - dep*0.6;
    
    % time vector setup
    start = 700;
    ending = 1800;
    maxt = ending-start;
    maxt = L;
    T = ones(length(dep),1)*(0:maxt);
    p = [c d]; % [-0.025, 0.2];
    
    S = twave_generator_full(T, dep, rep, p);
    tmp = S(1,:);
    sig_sim = A64 * S;
end

function [L12, idx, names] = get_lead_matrix()

    idx.V1 = 1;
    idx.V2 = 2;
    idx.V3 = 3;
    idx.V4 = 4;
    idx.V5 = 5;
    idx.V6 = 6;
    idx.RA = 7;
    idx.LA = 8;
    idx.LL = 9;

    L12 = zeros(12,64);

    % Lead I
    L12(1,idx.LA) =  1;
    L12(1,idx.RA) = -1;
    
    % Lead II
    L12(2,idx.LL) =  1;
    L12(2,idx.RA) = -1;
    
    % Lead III
    L12(3,idx.LL) =  1;
    L12(3,idx.LA) = -1;
    
    % aVR
    L12(4,idx.RA) =  1;
    L12(4,idx.LA) = -0.5;
    L12(4,idx.LL) = -0.5;
    
    % aVL
    L12(5,idx.LA) =  1;
    L12(5,idx.RA) = -0.5;
    L12(5,idx.LL) = -0.5;
    
    % aVF
    L12(6,idx.LL) =  1;
    L12(6,idx.RA) = -0.5;
    L12(6,idx.LA) = -0.5;
    
    % V1–V6
    for k = 1:6
        L12(6+k, idx.(['V' num2str(k)])) = 1;
        L12(6+k, idx.RA) = -1/3;
        L12(6+k, idx.LA) = -1/3;
        L12(6+k, idx.LL) = -1/3;
    end

    % names
    names = ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"];

end

function [dep, rep] = getTimesFromS(S)
    dep = first_crossing_thr(S, 0.5, true);
    rep = first_crossing_thr(S, 0.5, false);
end