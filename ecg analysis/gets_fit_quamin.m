close all

ecg = loadmat("data/VMedians/Normal.ecg_1.vmedianecg");
ecg_uni = ecg(7:end,:);

twave_start = 498;
twave_end = 734; % 720
twave = ecg_uni(:,twave_start:twave_end);

TT = rms(twave);
downslope_measured = 1 - cumsum(TT)/max(cumsum(TT));
[~, y_measLength] = size(downslope_measured);
y_meas = zeros(1, 500);
y_meas(1:y_measLength) = downslope_measured;

[theta_best, y_best, err, iter] = twave_optimize_quamin(y_meas);
fprintf("Best theta = [%.3f, %.3f, %.3f] \nError = %f\nIterations proceed = %d\n", theta_best(1), theta_best(2), theta_best(3), err, iter);

figure;
plot(y_meas(:), 'k', 'LineWidth', 1.2); hold on;
plot(y_best(:), 'r', 'LineWidth', 1.2);
legend( ...
    'real signal', ...
    sprintf("fitted -> theta=[%.3f, %.3f, %.3f] err=%.2f", theta_best(1), theta_best(2), theta_best(3), err) ...
    );
title('getS parameters Marquardt fit')
grid("on")


%%

% 1. Załadowanie całego datasetu
dataset_path = 'data/Dataset/';
loader = DatasetLoader(dataset_path);
loader.load_data();

% 2. Zdefiniowanie problemu optymalizacyjnego
initial_guess = [0.5, 0.5, 0.5]; % 3 parametry do Twojej nowej funkcji
modality = '12-lead';              % Modalność docelowa ('12-lead' lub '64-lead')

% 3. Utworzenie obiektu optymalizatora i przekazanie Twojej funkcji
optimizer = QuaminOptimizer(@get_y, initial_guess);

% 4. Pętla po wszystkich pacjentach
num_patients = loader.get_num_patients();

for i = 1:num_patients
    patient = loader.get_patient(i);

    
    % Pomijanie pacjentów z brakującymi danymi dla wybranej modalności
    [A, y_ref] = patient.get_modality_data(modality);
    if isempty(A) || isempty(y_ref)
        fprintf('Skiping patient %s (no data %s)\n', patient.PatientID, modality);
        continue; 
    end
    
    % Uruchomienie procesu optymalizacji dla konkretnego pacjenta
    fprintf('patient: %s...\n', patient.PatientID);
    [theta_best, y_fit, err, iter] = optimizer.run(patient, modality);
    
    % Wyświetlenie wyników
    fprintf('Succes in %d iterations,  error: %.4f | Theta = [%.3f, %.3f, %.3f]\n', ...
            iter, err, theta_best(1), theta_best(2), theta_best(3));
            
    % Tutaj możesz np. zapisać y_fit do pliku albo wygenerować wykresy dla pacjenta...
end

% -------------- FUNKCJA Z DEFINICJĄ MODELU --------------
function sig_sim = get_y(p, dep, rep, L, A)    
    mode = getsMode.NoReptCorrection; 
    param = [0.995, 0, p(1), p(2), p(3)];
    T = ones(length(dep),1)*(1:L);
    S = gets(T, dep, rep, param, mode);
    sig_sim = A * S;
end