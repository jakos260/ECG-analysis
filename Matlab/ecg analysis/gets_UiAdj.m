clearvars

% patient = 'normal_male'; offset = 1;
patient = 'normal_young_male'; offset = 169;
[ventri_ver, ventri_tri] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/ventricle.tri'));
[thorax_ver, thorax_tri] = loadtri_ecgsim(append(DATA_PATH, 'ECGsim_data/', patient, '/model/thorax.tri'));
 
q = initQtripy();
uiAdj = tmpGeneratorUiAdjuster(@get_h, patient, offset, @(v) q.values(v));

[~, min_idx] = min(uiAdj.fig.UserData.dep);
focus = ventri_ver(min_idx, :);
sim_values = uiAdj.fig.UserData.t_sim;
ref_values = uiAdj.fig.UserData.rep;

% QTriplot init
q.reset();
q.disable_debounce();
q.set_panels_number(2, 1);

% reference repolarisation time
q.set_active_panel(2, 1);
q.marker(focus, 'red', 5);
q.surface(ventri_ver, ventri_tri);
q.transparency(0.3);
q.values(ref_values);
q.gradient_bins(10);
q.text("reference rep", [0.6, 0.85]);

% simulated repolarisation time
q.set_active_panel(1, 1);
q.marker(focus, 'black', 5);
q.surface(ventri_ver, ventri_tri);
q.transparency(0.3);
q.gradient_bins(10);
q.text("simulated rep", [0.1, 0.85]);

q.enable_debounce(0.02)
% q.background_color('white')


function [tmp, sig_sim] = get_h(p, dep, rep, L, A)    
    % mode = 5; param = [theta(1), 0, theta(2), 0, theta(3)];
    % mode = 6; param = [theta(1), 0, theta(2), 0, theta(3)];
    % mode = getsMode.Exp_Spline; param = [0.995, 0, p(1)*2*0.96, p(2)*2, p(3)-0.2];
    mode = getsMode.Two_dy_Spline; param = [0.995, 0, p(1), p(2), p(3)];

    T = ones(length(dep),1)*(1:L);
    S = gets(T, dep, rep, param, mode);
    tmp = S;
    sig_sim = A * S;
end