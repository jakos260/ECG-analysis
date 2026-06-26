% Parametry stymulacji i symulacji
HRs = [60, 90, 120];                  % Częstości rytmu serca (bpm)
cell_types = [1, 2, 3];               % 1=Epi, 2=Mid, 3=Endo
cell_names = {'Epicardial', 'Mid-myocardial', 'Endocardial'};
colors = {'#0072BD', '#D95319', '#EDB120'}; % Niebieski, Pomarańczowy, Żółty

num_beats = 10;   % Liczba pobudzeń do osiągnięcia stanu zbliżonego do ustalonego
ISO = 0;          % Brak izoprenaliny
Stim_amp = -52;   % Amplituda prądu stymulacji (A/F)
Stim_dur = 1;     % Czas trwania stymulacji (ms)

% Inicjalizacja okna wykresu
figure('Name', 'APD for different HR', 'Position', [100, 100, 800, 900]);

for i = 1:length(HRs)
    HR = HRs(i);
    CL = 60000 / HR; % Przeliczenie HR na długość cyklu w milisekundach
    
    % Generowanie wektorów stymulacji dla funkcji TenTusscher2_IKr
    Stim_I = repmat(Stim_amp, 1, num_beats);
    Stim_T = repmat(Stim_dur, 1, num_beats);
    Stim_Int = repmat(CL, 1, num_beats);
    
    % Utworzenie subplota dla danego HR
    subplot(3, 1, i);
    hold on;
    
    for j = 1:length(cell_types)
        CT = cell_types(j);
        
        % Wywołanie modelu
        [t, V] = TenTusscher2_IKr(Stim_I, Stim_T, ISO, Stim_Int, CT);
        
        % Wyznaczanie czasu rozpoczęcia ostatniego impulsu
        % Zgodnie z konstrukcją funkcji, impulsy następują po sobie co Stim_Int
        t_start_last_beat = sum(Stim_Int(1:end-1));
        
        % Wyciągamy dane tylko dla ostatniego pobudzenia
        idx = find(t >= t_start_last_beat & t <= t_start_last_beat + CL);
        
        % Przesunięcie wektora czasu, aby początek ostatniego impulsu wypadał w t=0
        t_plot = t(idx) - t_start_last_beat;
        V_plot = V(idx);
        
        % Rysowanie krzywej
        plot(t_plot, V_plot, 'Color', colors{j}, 'LineWidth', 1.5, 'DisplayName', cell_names{j});
    end
    
    % Formatowanie pojedynczego subplota
    title(sprintf('Heart Rate = %d bpm (Cycle Length = %.0f ms)', HR, CL), 'FontSize', 12);
    xlabel('time (ms)');
    ylabel('(mV)');
    
    % Stała oś X (np. od -20 do 450 ms) na wszystkich wykresach 
    % ułatwia bezpośrednie wizualne porównanie skracania APD
    xlim([-20, 500]); 
    ylim([-95, 45]);
    
    grid on;
    legend('Location', 'northeast');
end

% Dodanie wspólnego tytułu
sgtitle('AP (Ten Tusscher 2006) for different HR', 'FontSize', 14, 'FontWeight', 'bold');