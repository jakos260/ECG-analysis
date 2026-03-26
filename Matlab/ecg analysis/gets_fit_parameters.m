close all;

ecg = loadmat("data/VMedians/Normal.ecg_1.vmedianecg");
ecg_uni = ecg(7:end,:);
figure(1);
plot(rms(ecg), 'b');
hold on
plot(rms(ecg_uni), 'k');

twave_start = 498;
twave_end = 734; % 720
twave = ecg_uni(:,twave_start:twave_end);
plot(twave_start:twave_end, rms(twave), 'r', 'LineWidth', 2);
legend("ecg rms", "uni\_ecg rms", "t-wave");
grid on;

%% optimalization

% original signal
TT = rms(twave);
downslope_measured = 1 - cumsum(TT)/max(cumsum(TT));
[~, y_measLength] = size(downslope_measured);
y_meas = zeros(1, 500);
y_meas(1:y_measLength) = downslope_measured;

% signal generator
function y = generator(theta, t)
    p = [0.2, 0, theta(1), 0, theta(2)]; % theta(2)
    dep = 0;
    rep = theta(3);
    S = gets_v2(t, dep, rep, p, 4);
    sMaxIdx = 36;
    y = zeros(1, 500);
    y(1:500-sMaxIdx+1) = S(sMaxIdx:500); % length(y_meas) hardcoded !
end

function y = generator_5(theta, t)
    p = [0.2, 0, theta(1), 0, theta(2)]; % theta(2)
    dep = 0;
    rep = theta(3);
    S = gets_v2(t, dep, rep, p, 5);
    sMaxIdx = 36;
    y = zeros(1, 500);
    y(1:500-sMaxIdx+1) = S(sMaxIdx:500); % length(y_meas) hardcoded !
end

function S = new_generator(theta, t)

    % parameters
    a = theta(1);   % slope of downslope
    rep = theta(2); % repolarization time shift
    k  = theta(3);  % ST slope / angle

    % shift time axis
    ts = t - rep;

    % --- ST section ---
    L0 = 20;        % base ST length
    alpha = 120;    % scaling for angle
    L = L0 + alpha * k;

    ST = ones(size(t));
    idx = ts >= 0 & ts <= L;
    ST(idx) = 1 - (ts(idx) / L); 
    ST(ts > L) = 0;

    % --- Down-slope sigmoid ---
    down = 1 ./ (1 + exp(a * ts));

    % final waveform
    S = ST .* down;


end

function E = error_fun(theta, t, y_meas)
    y = generator(theta, t);

    % Dopasowanie długości
    N = min(length(y), length(y_meas));
    y = y(1:N);
    y_meas = y_meas(1:N);

    % Klasyczny błąd MSE
    E = sum((y - y_meas).^2);
end


function [theta_opt, y_fit, final_error] = optimize_T_wave(y_meas, t)

    % --- Zakres czasu sygnału — zakładam że y_meas jest długości 500 ---
    if nargin < 2
        t = 1:length(y_meas);
    end

    % --- Inicjalizacja parametrów ---
    % theta = [p3, p5, rep]
    theta0 = [-0.04, 0.05, 180];   % sensowne startowe
    lb     = [-0.1, 0.5, 100];
    ub     = [0, 1, 450];

    % --- Funkcja błędu ---
    costfun = @(theta) error_fun(theta, t, y_meas);

    % --- Optymalizacja ---
    options = optimset('Display','iter', 'MaxIter', 300);
    theta_opt = fmincon(costfun, theta0, [],[],[],[], lb, ub, [], options);

    % --- Rekonstrukcja sygnału dla optimum ---
    y_fit = generator_5(theta_opt, t);
    final_error = costfun(theta_opt);

end


t = (1:500);

theta0  = [-0.04, 0.9, 185];    % starting point [p3, p5, rep]
lb      = [-0.1, 0.5, 0];
ub      = [0, 1, 250];

theta0_5  = [-0.04, 0.05, 180];
lb_5      = [-0.1, 0.5, 100];
ub_5      = [0, 1, 450];

n_theta0 = [-0.05, 150, 0.7];
n_lb     = [-1,   0,   0];
n_ub     = [ 1, 300,   1];

% ___ lsq ___ -> good for continous functions
[theta_opt, resnorm, residual] = lsqcurvefit( ...
    @(theta, t) generator(theta, t), ...
    theta0, t, y_meas, lb, ub);

[n_theta_opt, n_resnorm, n_residual] = lsqcurvefit( ...
    @(n_theta, t) new_generator(n_theta, t), ...
    n_theta0, t, y_meas, n_lb, n_ub);

[theta_5_best, y_5_best, err_5] = optimize_T_wave(y_meas);

% ___ patternsearch ___ -> requires a GADS_Toolbox license
% theta_opt = optimoptions('patternsearch','Display','iter','UseCompletePoll',true);
% [theta,fval] = patternsearch(@(th) cost(th), theta0, [],[],[],[], lb, ub, [], theta_opt);


% fprintf("theta_opt = [%f, %f, %f]\n", theta_opt(1), theta_opt(2), theta_opt(3));
% fprintf("resnorm = %f", resnorm);


%%

downslope_4 = generator(theta_opt, t);
new_downslope = new_generator(n_theta_opt, t);

figure(2);
plot(downslope_measured, 'r');
hold on
plot(new_downslope, 'k');
hold on
plot(downslope_4, 'b');
hold on
plot(y_5_best, 'g');
hold on

grid on
legend( ...
    "measured", ...
    sprintf("cn gen     theta = (%.3f, %.3f, %.3f) err=%f", n_theta_opt(1), n_theta_opt(2), n_theta_opt(3), n_resnorm), ...
    sprintf("my gen (4) theta = (%.3f, %.3f, %.3f) err=%f", theta_opt(1), theta_opt(2), theta_opt(3), resnorm), ...
    sprintf("my gen (5) theta = (%.3f, %.3f, %.3f) err=%f", theta_5_best(1), theta_5_best(2), theta_5_best(3), err_5) ...
);

figure(3);

plot(diff(downslope_measured), 'r');
hold on
plot(diff(new_downslope), 'k');
hold on
plot(diff(downslope_4), 'b');
hold on
plot(diff(y_5_best), 'g');
hold on
% plot(diff(downslope_sim), 'k');
grid on
legend( ...
    "measured", ...
    sprintf("cn gen     err=%f", n_resnorm), ...
    sprintf("my gen (4) err=%f", resnorm), ...
    sprintf("my gen (5) err=%f", resnorm_5) ...
);

    % sprintf("my theta = (%.3f, %.3f, %d)", my_theta(1), my_theta(2), my_theta(3)), ...



