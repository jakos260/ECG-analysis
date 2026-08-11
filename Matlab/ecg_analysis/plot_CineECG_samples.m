path = 'C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data\Dataset\';

close all

L = 2; % 1:99 lead number
for i = 1:10
    subject = sprintf('IKEM_Pat%03d', i);
    figure('Name', subject, 'NumberTitle', 'off');
    EcgDataPath = fullfile(path, subject, 'signals', 'ECG_DATA');
    CineEcgResultPath = fullfile(path, subject, 'signals', sprintf('IKEM_Pat%03d.iecg', i));
    cine_ecg = readCineEcgResults(CineEcgResultPath,'ecgdir', EcgDataPath, 'domedian');
    for n=1:size(cine_ecg.MEDIANDATA.VENTRICULAR, 2)
        signal  = cine_ecg.MEDIANDATA.VENTRICULAR{n}.beats.ECGbeat;
        start_t = cine_ecg.MEDIANDATA.VENTRICULAR{n}.beats.onsetqrs;
        cut_t   = cine_ecg.MEDIANDATA.VENTRICULAR{n}.beats.qrstduration;
        end_t   = min(size(signal, 2), start_t+cut_t);
    
        % diff_2 = abs(diff(signal(:,:), 2, 2));
        % stim_signal = max(diff_2);
        % [stim_signal_max_y, stim_signal_max_x] = max(stim_signal);
        % 
        % start_t = stim_signal_max_x+3;
        % cut_t = size(signal, 2) - start_t;
        
        subplot(3,3,n)
        plot(signal(L,:), 'k');
        hold on
        t = linspace(start_t, start_t+cut_t, cut_t+1);
        plot(baselinecor(cine_ecg.MEDIANDATA.VENTRICULAR{n}.ECG(L,1:1000))+0.5, 'b');
        title(cine_ecg.MEDIANDATA.VENTRICULAR{n}.ECGfilename, 'Interpreter', 'none');
        if end_t == start_t+cut_t
            plot(t, signal(L, start_t:end_t), 'r');
        else
            plot(start_t+cut_t, 0.1, 'r*');
            plot(start_t, 0.1, 'r*');
        end
    end
end