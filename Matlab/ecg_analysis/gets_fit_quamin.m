dataset_path = 'data/Dataset/';
loader = DatasetLoader(dataset_path);
loader.load_data();

initial_guess = [0.1, 0.6, 0.2]; % Start parameters for get_y function
LB = [ 0.0, 0.0, 0.0 ]; % Lower Bounds
UB = [ 1.0, 1.0, 1.0 ]; % Upper Bounds
modality = '12-lead';

optimizer = QuaminOptimizer(@get_y, initial_guess, LB, UB);
[theta_best, global_error, iterations] = optimizer.run(loader, modality);

fprintf('Best Theta = [%.3f, %.3f, %.3f]\n', theta_best(1), theta_best(2), theta_best(3));

patient = loader.Patients(1);
y_best = rms(get_y(theta_best, patient.dep, patient.rep, patient.L, patient.A12));
figure;
plot(rms(patient.y12), 'k', 'LineWidth', 1.2); hold on;
plot(y_best, 'r', 'LineWidth', 1.2);
legend( ...
    'real signal', ...
    sprintf("Best Theta=[%.3f, %.3f, %.3f]", theta_best(1), theta_best(2), theta_best(3)) ...
    );
title('getS parameters')
grid("on")

% -------------- MODEL DEFINITION --------------
function sig_sim = get_y(p, dep, rep, L, A)    
    mode = getsMode.Two_dy_Spline; 
    param = [0.995, 0, p(1), p(2), p(3)];
    T = ones(length(dep),1)*(1:L);
    S = gets(T, dep, rep, param, mode);
    sig_sim = A * S;
end
